"""Audit of the phase in a VCF, and cis/trans scenarios.

GEMPROT accepts a phase only when the data shows it. Two `|` heterozygotes are phased with each
other only if they are in the same `PS` block. If the VCF has no PS field, the phase is assumed
global, as with statistical or trio phasing. An unphased `/` genotype or a separate block can be
flipped. Each possible combination is a **scenario**.

haplo must not use a phase that the data does not show. Thus each copy of each scenario goes into
an intermediate VCF as a homozygous pseudo-sample (`a|a`). Then each protein maps to one known copy.
"""

from __future__ import annotations

import gzip
import itertools
from dataclasses import dataclass, field
from pathlib import Path

MAX_SCENARIOS = 16


# --- VCF reading --------------------------------------------------------------------------


@dataclass(frozen=True)
class Genotype:
    alleles: tuple[int | None, ...]  # None = missing allele. 1 item if haploid
    phased: bool
    ps: str | None  # phase block. None if absent or '.'

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
    meta: list[str]  # ## lines (FORMAT lines removed, GEMPROT writes them again)
    samples: list[str]
    records: list[Record]
    has_ps: bool
    ps_desc: dict[str, str] = field(default_factory=dict)  # PS block → description (fused phase)


def _parse_gt(gt: str, ps: str | None) -> Genotype:
    phased = "|" in gt
    sep = "|" if phased else "/"
    alleles = tuple(None if a == "." else int(a) for a in gt.split(sep))
    if len(alleles) == 1:
        phased = True  # haploid: nothing to phase
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


# --- audit per sample -----------------------------------------------------------------------


@dataclass(frozen=True)
class Call:
    """An alternative allele that the sample carries, and the copies that carry it (base orientation)."""

    record: int  # index in Vcf.records
    allele: int
    copies: frozenset[int]  # 0 = left of '|', 1 = right
    unit: int | None  # phase unit (None: homozygous or haploid, nothing to orient)
    label: str


@dataclass(frozen=True)
class PairPhase:
    a: Call
    b: Call
    relation: str  # cis, trans, homozygous, unknown
    reason: str


@dataclass
class SampleAudit:
    sample: str
    ploidy: int
    calls: list[Call]
    units: list[str]  # description of each phase unit
    pairs: list[PairPhase]
    notes: list[str] = field(default_factory=list)
    # base genotype per variant: alleles of copy 0 and copy 1
    base: dict[int, tuple[int, ...]] = field(default_factory=dict)
    unit_of: dict[int, int] = field(default_factory=dict)
    # missing alleles: (variant, copy in the VCF orientation, tested allele, name)
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
            what = "missing" if all(a is None for a in gt.alleles) else "partly missing"
            if len(rec.alts) == 1:
                notes.append(f"{rec.label(1)}: genotype {what}. GEMPROT tests presence and absence")
                missing += [(i, c, 1, rec.label(1)) for c, a in enumerate(gt.alleles) if a is None]
            else:
                notes.append(
                    f"{rec.label(1)}: genotype {what} at a multi-allelic site. "
                    "GEMPROT assumes the reference allele on the unknown copy"
                )
        if not gt.carried and None not in gt.alleles:
            continue
        alleles = tuple(a or 0 for a in gt.alleles)
        base[i] = alleles
        het = len(set(gt.alleles)) > 1
        unit = None
        if het:
            if not gt.phased:
                key, desc = ("unphased", i), f"{rec.label(max(gt.carried))} (unphased)"
            elif vcf.has_ps and gt.ps is None:
                key, desc = ("unphased", i), f"{rec.label(max(gt.carried))} (phased, no PS block)"
            elif vcf.has_ps:
                key, desc = ("ps", gt.ps), vcf.ps_desc.get(gt.ps, f"PS block {gt.ps}")
            else:
                key, desc = ("global",), "VCF phase (no PS field: assumed global)"
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
            pairs.append(PairPhase(a, b, "cis", "haploid"))
        elif len(a.copies) == ploidy or len(b.copies) == ploidy:
            pairs.append(PairPhase(a, b, "homozygous", "one variant is on both copies"))
        elif a.unit == b.unit:
            rel = "cis" if a.copies & b.copies else "trans"
            pairs.append(PairPhase(a, b, rel, unit_desc[a.unit]))
        else:
            pairs.append(PairPhase(a, b, "unknown", f"{unit_desc[a.unit]} / {unit_desc[b.unit]}"))
    return SampleAudit(sample, ploidy, calls, unit_desc, pairs, notes, base, unit_of, missing)


# --- scenarios ------------------------------------------------------------------------------


@dataclass(frozen=True)
class Scenario:
    index: int
    flipped: frozenset[int]  # units flipped relative to the VCF orientation
    copies: tuple[dict[int, int], ...]  # per copy: variant index → allele
    assumptions: tuple[str, ...] = ()  # assumptions about missing alleles

    def relation(self, a: Call, b: Call) -> str:
        ca = {c for c, alleles in enumerate(self.copies) if alleles.get(a.record) == a.allele}
        cb = {c for c, alleles in enumerate(self.copies) if alleles.get(b.record) == b.allele}
        return "cis" if ca & cb else "trans"


def scenarios(audit: SampleAudit, max_scenarios: int = MAX_SCENARIOS) -> list[Scenario]:
    if audit.n_scenarios > max_scenarios:
        raise ValueError(
            f"{audit.sample}: {len(audit.units)} independent phase units, "
            f"{audit.n_scenarios} scenarios (> {max_scenarios}). Phase this gene first"
        )
    out = []
    others = range(1, len(audit.units))  # the 1st unit sets the orientation (symmetry)
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
            assumptions.append(f"{label} {'present' if present else 'absent'} on copy {final + 1}")
        copies: tuple[dict[int, int], ...] = tuple({} for _ in range(audit.ploidy))
        for rec, alleles in base.items():
            if audit.unit_of.get(rec) in flipped:
                alleles = alleles[::-1]
            for c, allele in enumerate(alleles):
                copies[c][rec] = allele
        out.append(Scenario(k, flipped, copies, tuple(assumptions)))
    return out


# --- pseudo-sample VCF for haplo --------------------------------------------------------------

SEP = "::"


def pseudo_name(sample: str, scenario: int, copy: int) -> str:
    return f"{sample}{SEP}s{scenario}{SEP}c{copy}"


def parse_pseudo_name(name: str) -> tuple[str, int, int]:
    sample, s, c = name.rsplit(SEP, 2)
    return sample, int(s[1:]), int(c[1:])


def write_scenario_vcf(vcf: Vcf, plan: dict[str, list[Scenario]], out: Path) -> Path:
    """One column per (sample, scenario, copy), with a homozygous `a|a` genotype."""
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
