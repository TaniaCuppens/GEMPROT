"""Variants retenus pour les scénarios de phase : ceux qui peuvent changer la protéine.

Un variant hors CDS (intron, UTR, amont) ne change pas la séquence traduite par haplo : hétérozygote
non phasé, il double le nombre de scénarios sans changer leurs protéines. Un synonyme non plus,
sauf s'il partage un codon avec un autre variant (effet combiné).

On garde donc les sites dont un allèle est codant et change la protéine, plus les synonymes qui
partagent un codon avec un autre site codant. Les autres sont listés à part ; les variants
d'épissage y sont signalés (effet non modélisé).

Le tri a lieu après la fusion de la phase : un variant intronique phasé peut relier deux blocs.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from .annotate import IsolatedEffect, VepResult, sites_id
from .phase import Record, Vcf
from .phasing import FusedPhase

SPLICE = ("splice_acceptor_variant", "splice_donor_variant", "splice_donor_5th_base_variant")


@dataclass(frozen=True)
class Outside:
    """Site écarté des scénarios."""

    record: Record
    consequence: str  # conséquence VEP sur le transcrit, ou « hors transcrit »

    @property
    def splice(self) -> bool:
        return self.consequence in SPLICE


def _effects(vcf: Vcf, vep: VepResult, index: int) -> list[IsolatedEffect]:
    alleles = range(1, len(vcf.records[index].alts) + 1)
    return [vep.effects[key] for allele in alleles if (key := sites_id((index, allele))) in vep.effects]


def _neutral(effect: IsolatedEffect) -> bool:
    return effect.hgvsp is not None and effect.hgvsp.endswith("=")


def protein_records(vcf: Vcf, vep: VepResult) -> tuple[list[int], list[Outside]]:
    """Indices des sites gardés, et sites écartés."""
    effects = {index: _effects(vcf, vep, index) for index in range(len(vcf.records))}
    coding = {
        index: [effect for effect in site_effects if effect.protein_start is not None]
        for index, site_effects in effects.items()
    }
    sites_at_codon: dict[int, set[int]] = {}
    for index, site_effects in coding.items():
        for effect in site_effects:
            sites_at_codon.setdefault(effect.protein_start, set()).add(index)

    def changes_protein(index: int) -> bool:
        return any(not _neutral(effect) or len(sites_at_codon[effect.protein_start]) > 1 for effect in coding[index])

    keep = [index for index in effects if changes_protein(index)]
    outside = [
        Outside(vcf.records[index], effects[index][0].consequence if effects[index] else "hors transcrit")
        for index in effects
        if not changes_protein(index)
    ]
    return keep, outside


def subset(vcf: Vcf, keep: list[int]) -> Vcf:
    return Vcf(vcf.meta, vcf.samples, [vcf.records[index] for index in keep], vcf.has_ps, vcf.ps_desc)


def _new_indices(keep: list[int]) -> dict[int, int]:
    return {old: new for new, old in enumerate(keep)}


def remap_vep(vep: VepResult, keep: list[int]) -> VepResult:
    new_index = _new_indices(keep)
    parsed = {site_id: tuple(int(part) for part in site_id[1:].split("_")) for site_id in vep.effects}
    effects = {
        sites_id((new_index[index], allele)): vep.effects[site_id]
        for site_id, (index, allele) in parsed.items()
        if index in new_index
    }
    return VepResult(effects, vep.protein_id, vep.mane_ensembl)


def remap_phase(fused: FusedPhase, keep: list[int]) -> FusedPhase:
    new_index = _new_indices(keep)
    pairs = [
        replace(pair, a=new_index[pair.a], b=new_index[pair.b])
        for pair in fused.pairs
        if pair.a in new_index and pair.b in new_index
    ]
    orientation = {new_index[index]: value for index, value in fused.orientation.items() if index in new_index}
    return replace(fused, pairs=pairs, orientation=orientation)
