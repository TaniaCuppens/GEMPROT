"""Variant selection before the scenarios: only the variants that change the protein count."""

from pathlib import Path

import pytest

from gemprot import annotate as an
from gemprot import phasing as ph
from gemprot import scope
from gemprot.phase import Genotype, Record, Vcf, audit_sample, read_vcf

DATA = Path(__file__).parent / "data"
I114T, R197Q, R268K = 18400344, 18400593, 18400806


@pytest.fixture(scope="module")
def nat2():
    vcf = read_vcf(DATA / "nat2_trio.grch38.vcf")
    return vcf, an.parse_vep(DATA / "annot" / "nat2_trio.vep.json", "NM_000015.3")


def test_nat2_keeps_missense_only(nat2):
    vcf, vep = nat2
    keep, outside = scope.protein_records(vcf, vep)
    assert [vcf.records[index].pos for index in keep] == [I114T, R197Q, R268K]
    consequences = {site.record.pos: site.consequence for site in outside}
    assert consequences[18400484] == "synonymous_variant" and consequences[18399513] == "intron_variant"
    assert not any(site.splice for site in outside)


def test_nat2_scenarios_drop(nat2):
    vcf, vep = nat2
    keep, _ = scope.protein_records(vcf, vep)
    kept = scope.subset(vcf, keep)
    assert audit_sample(vcf, "NA12891").n_scenarios == 256
    assert audit_sample(kept, "NA12891").n_scenarios == 4
    assert audit_sample(kept, "NA12878").n_scenarios == 2


def test_phase_fused_before_selection_is_kept(nat2):
    vcf, vep = nat2
    family = ph.read_ped(DATA / "nat2_trio.ped")
    trio, _ = ph.trio_orientations(vcf, "NA12878", *family["NA12878"])
    fused = ph.fuse(vcf, "NA12878", {"trio": trio})
    keep, _ = scope.protein_records(vcf, vep)
    kept = scope.subset(ph.apply(vcf, {"NA12878": fused}), keep)
    remapped = scope.remap_phase(fused, keep)
    assert audit_sample(kept, "NA12878").n_scenarios == 1
    (pair,) = remapped.pairs
    assert (kept.records[pair.a].pos, kept.records[pair.b].pos) == (I114T, R268K)
    assert pair.final == "trans"


def test_remap_vep(nat2):
    vcf, vep = nat2
    keep, _ = scope.protein_records(vcf, vep)
    remapped = scope.remap_vep(vep, keep)
    assert remapped.effects[an.sites_id((2, 1))].hgvsp == "p.Arg268Lys"
    assert len(remapped.effects) == 3


# --- built cases ------------------------------------------------------------------------------


def _vcf(positions):
    genotype = Genotype((0, 1), False, None)
    records = [Record("1", pos, f"v{pos}", "A", ["G"], {"S": genotype}) for pos in positions]
    return Vcf([], ["S"], records, False)


def _effect(consequence, hgvsp=None, codon=None):
    return an.IsolatedEffect(None, hgvsp, consequence, codon)


def test_synonymous_sharing_codon_is_kept():
    vcf = _vcf([100, 101, 200])
    vep = an.VepResult({
        "v0_1": _effect("missense_variant", "p.Lys10Arg", 10),
        "v1_1": _effect("synonymous_variant", "p.Lys10=", 10),
        "v2_1": _effect("synonymous_variant", "p.Leu40=", 40),
    })
    keep, outside = scope.protein_records(vcf, vep)
    assert keep == [0, 1] and [site.record.pos for site in outside] == [200]


def test_splice_variant_is_flagged():
    vcf = _vcf([100, 300])
    vep = an.VepResult({"v0_1": _effect("missense_variant", "p.Lys10Arg", 10),
                        "v1_1": _effect("splice_donor_variant")})
    _, (site,) = scope.protein_records(vcf, vep)
    assert site.splice and site.consequence == "splice_donor_variant"


def test_variant_outside_transcript():
    vcf = _vcf([100, 900])
    vep = an.VepResult({"v0_1": _effect("missense_variant", "p.Lys10Arg", 10)})
    _, (site,) = scope.protein_records(vcf, vep)
    assert site.consequence == "outside transcript"
