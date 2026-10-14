"""Sources de phase et fusion des preuves.

Chaque source dit, pour des paires d'hétérozygotes d'un individu, si les allèles sont en cis ou en trans :

- **trio** : transmission mendélienne, à partir des génotypes des parents présents dans le VCF et
  d'un fichier `.ped` ; calculée ici, sans outil ;
- **reads** : WhatsHap sur les BAM/CRAM (blocs `PS`) ;
- **vcf** : la phase fournie en entrée (origine inconnue : statistique, trio ou reads).

Les sources restent séparées jusqu'à la fusion pour que leurs désaccords restent visibles. Fusion :
par ordre de fiabilité (trio > reads > vcf), chaque relation est ajoutée à un graphe (union-find avec
parité) ; une relation contredite par une source plus fiable est rejetée et signalée. Une source peut
ainsi relier des blocs d'une autre (le trio joint deux blocs de reads). Le résultat est réécrit en
GT/PS, un bloc par composante, puis passe par l'audit habituel (`phase.audit_sample`).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Sequence

from .config import Config
from .phase import Genotype, Record, Vcf, read_vcf
from .tools import run_tool

SOURCES = ("trio", "reads", "vcf")  # ordre de fiabilité
SOURCE_LABELS = {"trio": "trio", "reads": "reads (WhatsHap)", "vcf": "VCF d'entrée"}


@dataclass(frozen=True)
class Orientation:
    """Ce qu'une source dit d'un site hétérozygote : son bloc, et si la copie 0 porte le plus grand allèle."""

    block: str
    flip: bool


def het_sites(vcf: Vcf, sample: str) -> dict[int, tuple[int, int]]:
    """Sites hétérozygotes diploïdes complets → (petit allèle, grand allèle)."""
    out = {}
    for i, rec in enumerate(vcf.records):
        a = rec.genotypes[sample].alleles
        if len(a) == 2 and None not in a and a[0] != a[1]:
            out[i] = (min(a), max(a))
    return out


def orientations_from_vcf(vcf: Vcf, sample: str) -> dict[int, Orientation]:
    out = {}
    for i in het_sites(vcf, sample):
        gt = vcf.records[i].genotypes[sample]
        if not gt.phased or (vcf.has_ps and gt.ps is None):
            continue
        out[i] = Orientation(gt.ps if vcf.has_ps else "global", gt.alleles[0] > gt.alleles[1])
    return out


# --- trio ------------------------------------------------------------------------------------


def read_ped(path: Path) -> dict[str, tuple[str | None, str | None]]:
    """Fichier .ped (famille, individu, père, mère, …) → individu : (père, mère) ; `0` = inconnu."""
    out = {}
    for line in Path(path).read_text().splitlines():
        f = line.split()
        if len(f) < 4 or line.startswith("#"):
            continue
        out[f[1]] = (None if f[2] == "0" else f[2], None if f[3] == "0" else f[3])
    return out


def trio_orientations(
    vcf: Vcf, child: str, father: str | None, mother: str | None
) -> tuple[dict[int, Orientation], list[str]]:
    """Copie 0 = allèle transmis par le père. Un site n'est orienté que si une seule transmission est possible."""
    out, notes = {}, []
    for i, (lo, hi) in het_sites(vcf, child).items():
        rec = vcf.records[i]

        def possible(parent: str | None) -> set[int] | None:
            if parent is None or parent not in rec.genotypes:
                return None
            alleles = rec.genotypes[parent].alleles
            return None if None in alleles else set(alleles)

        pat, mat = possible(father), possible(mother)
        options = [(p, m) for p, m in ((lo, hi), (hi, lo))
                   if (pat is None or p in pat) and (mat is None or m in mat)]
        if not options:
            notes.append(f"{rec.label(hi)} : incompatibilité mendélienne avec les parents (erreur de génotype ou de novo)")
        elif len(options) == 1:
            out[i] = Orientation("trio", options[0][0] > options[0][1])
    return out, notes


# --- reads (WhatsHap) --------------------------------------------------------------------------


def bam_samples(bam: Path, cfg: Config) -> list[str]:
    header = run_tool(["samtools", "view", "-H", bam], cfg.vep.container, (*cfg.binds, bam.parent)).stdout
    names = []
    for line in header.splitlines():
        if line.startswith("@RG"):
            for tag in line.split("\t"):
                if tag.startswith("SM:") and tag[3:] not in names:
                    names.append(tag[3:])
    return names


def run_whatshap(vcf: Path, bams: Sequence[Path], samples: Sequence[str], out: Path, cfg: Config) -> Path:
    """Phasage par reads seuls (pas de --ped : le trio est une source à part, comparée ensuite)."""
    cmd = ["whatshap", "phase", "-o", out]
    if cfg.fasta:
        cmd += ["--reference", cfg.fasta]
    for s in samples:
        cmd += ["--sample", s]
    cmd += [vcf, *bams]
    binds = (*cfg.binds, vcf.parent, out.parent, *(b.parent for b in bams))
    if cfg.fasta:
        binds += (cfg.fasta.parent,)
    run_tool(cmd, cfg.vep.container, binds)
    return out


def _site_key(rec: Record) -> tuple:
    return rec.chrom, rec.pos, rec.ref, tuple(rec.alts)


def remap(other: Vcf, vcf: Vcf, sample: str) -> dict[int, Orientation]:
    """Orientations d'un autre VCF (sortie WhatsHap), ramenées aux indices de `vcf`."""
    index = {_site_key(r): i for i, r in enumerate(vcf.records)}
    het = het_sites(vcf, sample)
    out = {}
    for j, o in orientations_from_vcf(other, sample).items():
        i = index.get(_site_key(other.records[j]))
        if i in het:
            out[i] = o
    return out


# --- fusion ------------------------------------------------------------------------------------


@dataclass
class PairEvidence:
    a: int  # indices des sites dans le VCF
    b: int
    label_a: str
    label_b: str
    by_source: dict[str, str]  # source → cis/trans (allèles alternatifs les plus grands de chaque site)
    final: str | None = None
    status: str = "inconnue"  # accord, désaccord, source unique, déduite, inconnue

    @property
    def retained(self) -> str | None:
        """Source dont la relation est retenue (la plus fiable qui se prononce)."""
        return next((s for s in SOURCES if s in self.by_source), None)


@dataclass
class FusedPhase:
    sample: str
    sources: list[str]  # sources disponibles pour cet individu
    pairs: list[PairEvidence]
    notes: list[str] = field(default_factory=list)
    orientation: dict[int, tuple[str, bool]] = field(default_factory=dict)  # site → (bloc, flip)
    block_desc: dict[str, str] = field(default_factory=dict)

    @property
    def disagreements(self) -> list[PairEvidence]:
        return [p for p in self.pairs if p.status == "désaccord"]


class _ParityUnionFind:
    def __init__(self) -> None:
        self.parent: dict[int, int] = {}
        self.parity: dict[int, bool] = {}  # parité par rapport au parent

    def find(self, x: int) -> tuple[int, bool]:
        if self.parent.setdefault(x, x) == x:
            self.parity.setdefault(x, False)
            return x, False
        root, p = self.find(self.parent[x])
        self.parent[x], self.parity[x] = root, self.parity[x] ^ p
        return root, self.parity[x]

    def union(self, a: int, b: int, parity: bool) -> bool:
        """Ajoute « flip(a) xor flip(b) = parity » ; False si contredit ce qui est déjà établi."""
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        if ra == rb:
            return (pa ^ pb) == parity
        self.parent[rb], self.parity[rb] = ra, pa ^ pb ^ parity
        return True


def fuse(vcf: Vcf, sample: str, evidence: dict[str, dict[int, Orientation]]) -> FusedPhase:
    het = het_sites(vcf, sample)
    sources = [s for s in SOURCES if evidence.get(s)]
    relations: dict[tuple[int, int], dict[str, str]] = {}
    for source in sources:
        by_block: dict[str, list[int]] = {}
        for i, o in evidence[source].items():
            by_block.setdefault(o.block, []).append(i)
        for sites in by_block.values():
            for a, b in itertools.combinations(sorted(sites), 2):
                same = evidence[source][a].flip == evidence[source][b].flip
                relations.setdefault((a, b), {})[source] = "cis" if same else "trans"

    uf = _ParityUnionFind()
    rejected = []
    for source in sources:
        for (a, b), rel in relations.items():
            if source in rel and not uf.union(a, b, rel[source] == "trans"):
                rejected.append((a, b, source))

    pairs = []
    for a, b in itertools.combinations(sorted(het), 2):
        rel = relations.get((a, b), {})
        (ra, pa), (rb, pb) = uf.find(a), uf.find(b)
        final = ("trans" if pa ^ pb else "cis") if ra == rb else None
        if len(set(rel.values())) > 1:
            status = "désaccord"
        elif len(rel) > 1:
            status = "accord"
        elif rel:
            status = "source unique"
        else:
            status = "déduite" if final else "inconnue"
        pairs.append(PairEvidence(a, b, vcf.records[a].label(het[a][1]), vcf.records[b].label(het[b][1]),
                                  rel, final, status))

    # désaccord indirect : la paire n'est contredite par aucune autre source, mais par une chaîne
    # de relations plus fiables (A-B et B-C donnent A-C)
    status_of = {(p.a, p.b): p.status for p in pairs}
    notes = [
        f"{pairs_label(vcf, het, a, b)} : {SOURCE_LABELS[source]} dit {relations[(a, b)][source]}, "
        f"contredit par d'autres relations plus fiables — non retenu"
        for a, b, source in rejected
        if status_of[(a, b)] != "désaccord"
    ]

    # blocs : composantes à au moins deux sites ; sources qui y ont contribué
    members: dict[int, list[int]] = {}
    for i in het:
        if i in uf.parent:
            members.setdefault(uf.find(i)[0], []).append(i)
    orientation, block_desc = {}, {}
    for k, (root, sites) in enumerate(sorted(members.items(), key=lambda kv: min(kv[1])), 1):
        if len(sites) < 2:
            continue
        block = f"F{k}"
        used = [s for s in sources if any(s in relations.get((a, b), {})
                                          for a, b in itertools.combinations(sorted(sites), 2))]
        block_desc[block] = f"bloc {k} ({' + '.join(SOURCE_LABELS[s] for s in used)})"
        for i in sites:
            orientation[i] = (block, uf.find(i)[1])
    return FusedPhase(sample, sources, pairs, notes, orientation, block_desc)


def pairs_label(vcf: Vcf, het: dict[int, tuple[int, int]], a: int, b: int) -> str:
    return f"{vcf.records[a].label(het[a][1])} / {vcf.records[b].label(het[b][1])}"


def apply(vcf: Vcf, fused: dict[str, FusedPhase]) -> Vcf:
    """Réécrit les hétérozygotes des individus fusionnés : `|` + bloc fusionné, sinon `/`."""
    records = []
    for i, rec in enumerate(vcf.records):
        gts = dict(rec.genotypes)
        for sample, f in fused.items():
            gt = gts[sample]
            a = gt.alleles
            if not (len(a) == 2 and None not in a and a[0] != a[1]):
                continue
            lo, hi = min(a), max(a)
            if i in f.orientation:
                block, flip = f.orientation[i]
                gts[sample] = Genotype((hi, lo) if flip else (lo, hi), True, block)
            else:
                gts[sample] = Genotype((lo, hi), False, None)
        records.append(replace(rec, genotypes=gts))
    desc = dict(vcf.ps_desc)
    for f in fused.values():
        desc.update(f.block_desc)
    if not vcf.has_ps:
        # individus non fusionnés : leur phase sans PS reste globale, comme dans le VCF d'origine
        desc["global"] = "phase du VCF (pas de champ PS : supposée globale)"
        for rec in records:
            for sample, gt in rec.genotypes.items():
                if sample not in fused and gt.phased and len(gt.alleles) > 1:
                    rec.genotypes[sample] = replace(gt, ps="global")
    return Vcf(vcf.meta, vcf.samples, records, True, desc)


def phase_sources(
    vcf: Vcf,
    vcf_path: Path,
    samples: Sequence[str],
    cfg: Config,
    workdir: Path,
    bams: Sequence[Path] = (),
    ped: Path | None = None,
) -> dict[str, FusedPhase]:
    evidence: dict[str, dict[str, dict[int, Orientation]]] = {s: {} for s in samples}
    notes: dict[str, list[str]] = {s: [] for s in samples}
    for s in samples:
        evidence[s]["vcf"] = orientations_from_vcf(vcf, s)
    if ped:
        family = read_ped(ped)
        for s in samples:
            if s in family and any(p in vcf.samples for p in family[s] if p):
                evidence[s]["trio"], notes[s] = trio_orientations(vcf, s, *family[s])
    if bams:
        with_reads = []
        for bam in bams:
            with_reads += [n for n in bam_samples(bam, cfg) if n in samples and n not in with_reads]
        if with_reads:
            out = read_vcf(run_whatshap(vcf_path, bams, with_reads, workdir / "whatshap.vcf", cfg))
            for s in with_reads:
                evidence[s]["reads"] = remap(out, vcf, s)
    fused = {}
    for s in samples:
        f = fuse(vcf, s, evidence[s])
        f.notes = notes[s] + f.notes
        fused[s] = f
    return fused
