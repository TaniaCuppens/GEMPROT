"""Effet de chaque haplotype sur la protéine, et verdict par individu.

On traduit la CDS de référence et celle de l'haplotype (fournie par haplo), puis on décrit la
différence en HGVS protéine. On repère aussi les effets combinés en cis :
- variant « masqué » : en aval d'un stop ou d'un décalage du cadre sur le même haplotype,
  il n'a plus d'effet propre (cas GJB2 R75Q + c.35delG, Iossa 2010) ;
- décalage du cadre compensé : deux indels en cis dont la somme garde le cadre.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .haplotypes import CdsDiff, CdsHaplotype, TranscriptHaplotypes

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AA[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa",
}


def translate(cds: str) -> str:
    """Traduit jusqu'au premier stop inclus ('*') ; sans stop, s'arrête au dernier codon complet."""
    prot = []
    for i in range(0, len(cds) - 2, 3):
        aa = CODON_TABLE.get(cds[i : i + 3].upper(), "X")
        prot.append(aa)
        if aa == "*":
            break
    return "".join(prot)


def _aa(seq: str, i: int) -> str:
    """Acide aminé en position 0-based i, en code à 3 lettres avec sa position 1-based."""
    return f"{THREE[seq[i]]}{i + 1}"


# --- différences en coordonnées CDS (notation haplo) ------------------------------------

_DIFF_RE = re.compile(
    r"^(?P<pos>\d+)(?:(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)|del(?P<del>[ACGTN]+|\{\d+\})|ins(?P<ins>[ACGTN]+|\{\d+\}))$"
)


def _seq_len(s: str) -> int:
    return int(s[1:-1]) if s.startswith("{") else len(s)


@dataclass(frozen=True)
class ParsedDiff:
    raw: CdsDiff
    pos: int  # position CDS 1-based
    length_change: int | None  # None si la notation n'est pas reconnue

    @property
    def codon(self) -> int:
        return (self.pos - 1) // 3 + 1

    @property
    def frameshifting(self) -> bool:
        return self.length_change is not None and self.length_change % 3 != 0


def parse_diff(d: CdsDiff) -> ParsedDiff:
    m = _DIFF_RE.match(d.diff)
    if not m:
        lead = re.match(r"\d+", d.diff)
        return ParsedDiff(d, int(lead.group()) if lead else 0, None)
    if m["ref"]:
        change = len(m["alt"]) - len(m["ref"])
    elif m["del"]:
        change = -_seq_len(m["del"])
    else:
        change = _seq_len(m["ins"])
    return ParsedDiff(d, int(m["pos"]), change)


def hgvs_c(d: CdsDiff, ref_cds: str) -> str:
    """Notation HGVS c. ; les délétions sont décalées en 3' comme le veut HGVS (haplo les aligne à gauche)."""
    m = _DIFF_RE.match(d.diff)
    if m and m["del"] and not m["del"].startswith("{"):
        start, seq = int(m["pos"]), m["del"]
        while start - 1 + len(seq) < len(ref_cds) and ref_cds[start - 1 + len(seq)] == seq[0]:
            start, seq = start + 1, seq[1:] + seq[0]
        end = start + len(seq) - 1
        return f"c.{start}del" if end == start else f"c.{start}_{end}del"
    if m and m["ref"]:
        return f"c.{m['pos']}{m['ref']}>{m['alt']}"
    return f"c.{d.diff}"


# --- effet protéique d'un haplotype -------------------------------------------------------


@dataclass(frozen=True)
class ProteinChange:
    kind: str  # missense, stop_gained, stop_lost, start_lost, inframe_indel, frameshift
    hgvs: str
    position: int  # premier acide aminé touché (1-based)


@dataclass
class HaplotypeEffect:
    haplotype: CdsHaplotype
    protein: str
    changes: list[ProteinChange] = field(default_factory=list)
    truncated_at: int | None = None  # dernier acide aminé normal avant la perte de la protéine
    masked: list[CdsDiff] = field(default_factory=list)  # variants sans effet propre (en aval)
    compensated_frameshift: bool = False

    @property
    def altered(self) -> bool:
        return bool(self.changes)

    @property
    def label(self) -> str:
        if not self.haplotype.diffs:
            return "référence"
        if not self.changes:
            return "synonyme"
        return ", ".join(c.hgvs for c in self.changes)


def _substitutions(ref: str, hap: str, offset: int = 0) -> list[ProteinChange]:
    return [
        ProteinChange("missense", f"p.{_aa(ref, offset + k)}{THREE[b]}", offset + k + 1)
        for k, (a, b) in enumerate(zip(ref, hap))
        if a != b
    ]


def _inframe(ref_p: str, hap_p: str, net: int = 0) -> list[ProteinChange]:
    """Différences sans décalage du cadre : substitutions, stop gagné/perdu, indels en phase.

    `net` : variation de longueur de la CDS (multiple de 3) ; non nulle → indel en phase.
    """
    if ref_p == hap_p:
        return []
    if hap_p[:1] != ref_p[:1] and ref_p[:1] == "M":
        return [ProteinChange("start_lost", "p.Met1?", 1)]
    k = len(ref_p) - 1
    if net == 0 and ref_p.endswith("*") and len(hap_p) > k and hap_p[:k] == ref_p[:k] and hap_p[k] != "*":
        ext = hap_p[k:]
        tail = f"*{len(ext) - 1}" if ext.endswith("*") else "*?"
        return [ProteinChange("stop_lost", f"p.Ter{k + 1}{THREE[ext[0]]}ext{tail}", k + 1)]
    # stop gagné : la protéine s'arrête plus tôt, le reste est comparé position à position
    if net == 0 and hap_p.endswith("*") and len(hap_p) < len(ref_p):
        k = len(hap_p) - 1
        changes = _substitutions(ref_p[:k], hap_p[:k])
        return changes + [ProteinChange("stop_gained", f"p.{_aa(ref_p, k)}*", k + 1)]
    if net == 0 and len(ref_p) == len(hap_p):
        return _substitutions(ref_p, hap_p)
    # indel en phase, ou stop perdu : on retire le préfixe et le suffixe communs
    i = 0
    while i < min(len(ref_p), len(hap_p)) and ref_p[i] == hap_p[i]:
        i += 1
    if ref_p[i:] == "*" and len(hap_p) > len(ref_p) - 1:
        ext = hap_p[i:]
        tail = f"*{len(ext) - 1}" if ext.endswith("*") else "*?"
        return [ProteinChange("stop_lost", f"p.Ter{i + 1}{THREE[ext[0]]}ext{tail}", i + 1)]
    s = 0
    while s < min(len(ref_p), len(hap_p)) - i and ref_p[-1 - s] == hap_p[-1 - s]:
        s += 1
    r, h = ref_p[i : len(ref_p) - s], hap_p[i : len(hap_p) - s]
    first = _aa(ref_p, i)
    span = first if len(r) <= 1 else f"{first}_{_aa(ref_p, i + len(r) - 1)}"
    if not h:
        return [ProteinChange("inframe_indel", f"p.{span}del", i + 1)]
    ins = "".join(THREE[a] for a in h)
    if not r:
        return [ProteinChange("inframe_indel", f"p.{_aa(ref_p, i - 1)}_{first}ins{ins}", i + 1)]
    return [ProteinChange("inframe_indel", f"p.{span}delins{ins}", i + 1)]


def haplotype_effect(hap: CdsHaplotype, ref_cds: str) -> HaplotypeEffect:
    ref_p, hap_p = translate(ref_cds), translate(hap.seq)
    effect = HaplotypeEffect(hap, hap_p)
    diffs = sorted((parse_diff(d) for d in hap.diffs), key=lambda p: p.pos)

    frameshifts = [d for d in diffs if d.frameshifting]
    net = len(hap.seq) - len(ref_cds)
    fs_codon = frameshifts[0].codon if frameshifts and net % 3 != 0 else None

    if fs_codon is None:
        effect.changes = _inframe(ref_p, hap_p, net)
        # deux indels en cis qui se compensent, sans stop entre les deux
        effect.compensated_frameshift = bool(frameshifts) and not any(
            c.kind == "stop_gained" for c in effect.changes
        )
        if effect.compensated_frameshift:
            # un seul événement entre les deux indels : décrit en delins (HGVS)
            lo, hi = frameshifts[0].codon, frameshifts[-1].codon
            inside = [c for c in effect.changes if lo <= c.position <= hi and c.kind == "missense"]
            if len(inside) > 1:
                a, z = inside[0].position, inside[-1].position
                ins = "".join(THREE[x] for x in hap_p[a - 1 : z])
                merged = ProteinChange(
                    "inframe_indel", f"p.{_aa(ref_p, a - 1)}_{_aa(ref_p, z - 1)}delins{ins}", a
                )
                effect.changes = sorted(
                    [c for c in effect.changes if c not in inside] + [merged], key=lambda c: c.position
                )
    else:
        # en amont du décalage : comparaison en phase ; à partir du décalage : un seul événement
        j = fs_codon - 1
        while j < min(len(ref_p), len(hap_p)) and ref_p[j] == hap_p[j]:
            j += 1  # le décalage peut ne changer l'acide aminé qu'un peu plus loin
        effect.changes = _inframe(ref_p[: fs_codon - 1], hap_p[: fs_codon - 1])
        if j < len(hap_p) and hap_p[j] == "*":
            effect.changes.append(ProteinChange("stop_gained", f"p.{_aa(ref_p, j)}*", j + 1))
        elif j < len(ref_p):
            tail = f"*{len(hap_p) - j}" if hap_p.endswith("*") else "*?"
            effect.changes.append(
                ProteinChange("frameshift", f"p.{_aa(ref_p, j)}{THREE[hap_p[j]]}fs{tail}", j + 1)
            )

    lost = [c.position for c in effect.changes if c.kind in ("stop_gained", "frameshift", "start_lost")]
    if lost:
        effect.truncated_at = min(lost) - 1
        cut = effect.truncated_at + 1  # 1er acide aminé perdu ; la cause est toujours en amont ou dessus
        effect.masked = [d.raw for d in diffs if d.codon > cut]
    return effect


# --- verdict par individu --------------------------------------------------------------------


@dataclass
class SampleVerdict:
    sample: str
    copies: list[HaplotypeEffect]

    @property
    def altered_copies(self) -> int:
        return sum(c.altered for c in self.copies)

    @property
    def summary(self) -> str:
        n, total = self.altered_copies, len(self.copies)
        if n == 0:
            return "aucune copie altérée"
        if n == total:
            return f"{'les deux copies' if total == 2 else 'la copie'} altérée{'s' if total == 2 else ''}"
        return f"{n} copie altérée sur {total} (une copie normale)"


def classify_pseudo_samples(tx: TranscriptHaplotypes) -> dict[str, dict[int, SampleVerdict]]:
    """Sortie haplo d'un VCF de pseudo-échantillons (phase.write_scenario_vcf) :
    individu → scénario → verdict, les copies dans l'ordre (gauche puis droite du '|')."""
    from .phase import parse_pseudo_name

    effects: dict[int, HaplotypeEffect] = {}
    copies: dict[tuple[str, int], dict[int, HaplotypeEffect]] = {}
    for name, haps in tx.sample_haplotypes().items():
        sample, scenario, copy = parse_pseudo_name(name)
        hap = haps[0]  # pseudo-échantillon homozygote : deux fois le même haplotype
        if id(hap) not in effects:
            effects[id(hap)] = haplotype_effect(hap, tx.ref_cds)
        copies.setdefault((sample, scenario), {})[copy] = effects[id(hap)]
    out: dict[str, dict[int, SampleVerdict]] = {}
    for (sample, scenario), by_copy in sorted(copies.items()):
        out.setdefault(sample, {})[scenario] = SampleVerdict(sample, [by_copy[c] for c in sorted(by_copy)])
    return out


def classify_transcript(tx: TranscriptHaplotypes) -> dict[str, SampleVerdict]:
    effects = {id(h): haplotype_effect(h, tx.ref_cds) for h in tx.haplotypes}
    return {
        sample: SampleVerdict(sample, [effects[id(h)] for h in haps])
        for sample, haps in sorted(tx.sample_haplotypes().items())
    }
