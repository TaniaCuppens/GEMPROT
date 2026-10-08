"""Audit de la phase d'un VCF et scénarios cis/trans.

Une phase n'est crue que si elle est démontrée : deux hétérozygotes `|` ne sont phasés entre eux
que s'ils sont dans le même bloc `PS` (ou si le VCF n'a pas de champ PS : phase supposée globale,
ce qui est le cas d'un phasage statistique ou par trio). Un génotype non phasé `/` ou un bloc à
part peut être retourné : chaque combinaison possible est un **scénario**.

Pour que haplo ne puisse pas réutiliser une phase non démontrée, chaque copie de chaque scénario
est écrite comme un pseudo-échantillon homozygote (`a|a`) dans un VCF intermédiaire ; on sait
ainsi exactement quelle protéine correspond à quelle copie.
"""

from __future__ import annotations

import gzip
import itertools
from dataclasses import dataclass, field
from pathlib import Path

MAX_SCENARIOS = 16


# --- lecture du VCF -----------------------------------------------------------------------


@dataclass(frozen=True)
class Genotype:
    alleles: tuple[int | None, ...]  # None = allèle manquant ; 1 élément si haploïde
    phased: bool
    ps: str | None  # bloc de phase ; None si absent ou '.'

    @property
    def carried(self) -> set[int]:
        return {a for a in self.alleles if a}


@dataclass
class Record:
    chrom: str
    pos: int
    id: str
    ref: str
    alts: list[str]
    genotypes: dict[str, Genotype]

    def label(self, allele: int) -> str:
        if self.id not in (".", "") and len(self.alts) == 1:
            return self.id
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alts[allele - 1]}"


@dataclass
class Vcf:
    meta: list[str]  # lignes ## (sans les FORMAT, réécrites)
    samples: list[str]
    records: list[Record]
    has_ps: bool


def _parse_gt(gt: str, ps: str | None) -> Genotype:
    phased = "|" in gt
    sep = "|" if phased else "/"
    alleles = tuple(None if a == "." else int(a) for a in gt.split(sep))
    if len(alleles) == 1:
        phased = True  # haploïde : rien à phaser
    return Genotype(alleles, phased, None if ps in (None, ".", "") else ps)


def read_vcf(path: Path | str) -> Vcf:
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    meta, samples, records, has_ps = [], [], [], False
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                if line.startswith("##FORMAT=<ID=PS,"):
                    has_ps = True
                if not line.startswith("##FORMAT="):
                    meta.append(line)
                continue
            if line.startswith("#"):
                samples = line.split("\t")[9:]
                continue
            f = line.split("\t")
            keys = f[8].split(":")
            gts = {}
            for sample, value in zip(samples, f[9:]):
                vals = dict(zip(keys, value.split(":")))
                gts[sample] = _parse_gt(vals.get("GT", "."), vals.get("PS"))
            records.append(Record(f[0], int(f[1]), f[2], f[3], f[4].split(","), gts))
    return Vcf(meta, samples, records, has_ps)


# --- audit par individu ---------------------------------------------------------------------


@dataclass(frozen=True)
class Call:
    """Un allèle alternatif porté par l'individu, et les copies qui le portent (orientation de base)."""

    record: int  # indice dans Vcf.records
    allele: int
    copies: frozenset[int]  # 0 = gauche du '|', 1 = droite
    unit: int | None  # unité de phase (None : homozygote ou haploïde, rien à orienter)
    label: str


@dataclass(frozen=True)
class PairPhase:
    a: Call
    b: Call
    relation: str  # cis, trans, homozygote, inconnue
    reason: str


@dataclass
class SampleAudit:
    sample: str
    ploidy: int
    calls: list[Call]
    units: list[str]  # description de chaque unité de phase
    pairs: list[PairPhase]
    notes: list[str] = field(default_factory=list)
    # génotype de base par variant : allèles de la copie 0 et de la copie 1
    base: dict[int, tuple[int, ...]] = field(default_factory=dict)
    unit_of: dict[int, int] = field(default_factory=dict)
    # allèles manquants : (variant, copie dans l'orientation du VCF, allèle testé, nom)
    missing: list[tuple[int, int, int, str]] = field(default_factory=list)

    @property
    def n_scenarios(self) -> int:
        return 2 ** (max(len(self.units) - 1, 0) + len(self.missing))


def audit_sample(vcf: Vcf, sample: str) -> SampleAudit:
    calls: list[Call] = []
    units: dict[tuple, int] = {}
    unit_desc: list[str] = []
    notes: list[str] = []
    base: dict[int, tuple[int, ...]] = {}
    unit_of: dict[int, int] = {}
    missing: list[tuple[int, int, int, str]] = []
    ploidy = 1

    for i, rec in enumerate(vcf.records):
        gt = rec.genotypes[sample]
        ploidy = max(ploidy, len(gt.alleles))
        if None in gt.alleles:
            what = "manquant" if all(a is None for a in gt.alleles) else "partiellement manquant"
            if len(rec.alts) == 1:
                notes.append(f"{rec.label(1)} : génotype {what}, présence et absence testées")
                missing += [(i, c, 1, rec.label(1)) for c, a in enumerate(gt.alleles) if a is None]
            else:
                notes.append(f"{rec.label(1)} : génotype {what} sur un site multi-allélique, copie inconnue supposée de référence")
        if not gt.carried and None not in gt.alleles:
            continue
        alleles = tuple(a or 0 for a in gt.alleles)
        base[i] = alleles
        het = len(set(gt.alleles)) > 1
        unit = None
        if het:
            if not gt.phased:
                key, desc = ("unphased", i), f"{rec.label(max(gt.carried))} (non phasé)"
            elif vcf.has_ps and gt.ps is None:
                key, desc = ("unphased", i), f"{rec.label(max(gt.carried))} (phasé sans bloc PS)"
            elif vcf.has_ps:
                key, desc = ("ps", gt.ps), f"bloc PS {gt.ps}"
            else:
                key, desc = ("global",), "phase du VCF (pas de champ PS : supposée globale)"
            if key not in units:
                units[key] = len(unit_desc)
                unit_desc.append(desc)
            unit = units[key]
            unit_of[i] = unit
        if not gt.carried:
            continue
        for allele in sorted(gt.carried):
            copies = frozenset(c for c, a in enumerate(alleles) if a == allele)
            calls.append(Call(i, allele, copies, unit, rec.label(allele)))

    pairs = []
    for a, b in itertools.combinations(calls, 2):
        if a.record == b.record:
            continue
        if ploidy == 1:
            pairs.append(PairPhase(a, b, "cis", "haploïde"))
        elif len(a.copies) == ploidy or len(b.copies) == ploidy:
            pairs.append(PairPhase(a, b, "homozygote", "un des variants est sur les deux copies"))
        elif a.unit == b.unit:
            rel = "cis" if a.copies & b.copies else "trans"
            pairs.append(PairPhase(a, b, rel, unit_desc[a.unit]))
        else:
            pairs.append(PairPhase(a, b, "inconnue", f"{unit_desc[a.unit]} / {unit_desc[b.unit]}"))
    return SampleAudit(sample, ploidy, calls, unit_desc, pairs, notes, base, unit_of, missing)


# --- scénarios ------------------------------------------------------------------------------


@dataclass(frozen=True)
class Scenario:
    index: int
    flipped: frozenset[int]  # unités retournées par rapport à l'orientation du VCF
    copies: tuple[dict[int, int], ...]  # par copie : indice de variant → allèle
    assumptions: tuple[str, ...] = ()  # hypothèses sur les allèles manquants

    def relation(self, a: Call, b: Call) -> str:
        ca = {c for c, alleles in enumerate(self.copies) if alleles.get(a.record) == a.allele}
        cb = {c for c, alleles in enumerate(self.copies) if alleles.get(b.record) == b.allele}
        return "cis" if ca & cb else "trans"


def scenarios(audit: SampleAudit, max_scenarios: int = MAX_SCENARIOS) -> list[Scenario]:
    if audit.n_scenarios > max_scenarios:
        raise ValueError(
            f"{audit.sample} : {len(audit.units)} unités de phase indépendantes, "
            f"{audit.n_scenarios} scénarios (> {max_scenarios}) ; il faut phaser ce gène"
        )
    out = []
    others = range(1, len(audit.units))  # la 1re unité fixe l'orientation (symétrie)
    flips = itertools.product((False, True), repeat=len(others))
    fills = list(itertools.product((False, True), repeat=len(audit.missing)))
    for k, (flip, fill) in enumerate(itertools.product(flips, fills)):
        flipped = frozenset(u for u, f in zip(others, flip) if f)
        base = {rec: list(alleles) for rec, alleles in audit.base.items()}
        assumptions = []
        for (rec, c, allele, label), present in zip(audit.missing, fill):
            if present:
                base[rec][c] = allele
            final = len(base[rec]) - 1 - c if audit.unit_of.get(rec) in flipped else c
            assumptions.append(f"{label} {'présent' if present else 'absent'} sur la copie {final + 1}")
        copies: tuple[dict[int, int], ...] = tuple({} for _ in range(audit.ploidy))
        for rec, alleles in base.items():
            if audit.unit_of.get(rec) in flipped:
                alleles = alleles[::-1]
            for c, allele in enumerate(alleles):
                copies[c][rec] = allele
        out.append(Scenario(k, flipped, copies, tuple(assumptions)))
    return out


# --- VCF de pseudo-échantillons pour haplo -----------------------------------------------------

SEP = "::"


def pseudo_name(sample: str, scenario: int, copy: int) -> str:
    return f"{sample}{SEP}s{scenario}{SEP}c{copy}"


def parse_pseudo_name(name: str) -> tuple[str, int, int]:
    sample, s, c = name.rsplit(SEP, 2)
    return sample, int(s[1:]), int(c[1:])


def write_scenario_vcf(vcf: Vcf, plan: dict[str, list[Scenario]], out: Path) -> Path:
    """Une colonne par (individu, scénario, copie), génotype homozygote `a|a`."""
    columns = [
        (pseudo_name(sample, sc.index, c), sc.copies[c])
        for sample, scs in plan.items()
        for sc in scs
        for c in range(len(sc.copies))
    ]
    with open(out, "w") as fh:
        for line in vcf.meta:
            fh.write(line + "\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                           + [name for name, _ in columns]) + "\n")
        for i, rec in enumerate(vcf.records):
            gts = [f"{alleles.get(i, 0)}|{alleles.get(i, 0)}" for _, alleles in columns]
            fh.write("\t".join([rec.chrom, str(rec.pos), rec.id, rec.ref, ",".join(rec.alts),
                                ".", ".", ".", "GT", *gts]) + "\n")
    return out
