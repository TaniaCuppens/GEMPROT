"""Sources de phase (trio, reads, VCF) et fusion, sur le trio CEPH 1463 (NAT2) et des cas construits."""

from pathlib import Path

import pytest

from gemprot import phasing as ph
from gemprot.phase import Genotype, Record, Vcf, audit_sample, read_vcf

DATA = Path(__file__).parent / "data"
CHILD = "NA12878"
I114T, L161, R268K = 18400344, 18400484, 18400806


@pytest.fixture(scope="module")
def nat2():
    return read_vcf(DATA / "nat2_trio.grch38.vcf")


def _idx(vcf, pos):
    return next(i for i, r in enumerate(vcf.records) if r.pos == pos)


def _pair(fused, vcf, p1, p2):
    a, b = sorted((_idx(vcf, p1), _idx(vcf, p2)))
    return next(p for p in fused.pairs if (p.a, p.b) == (a, b))


def _evidence(vcf, *, trio=True, reads=True, panel=False):
    ev = {}
    if trio:
        family = ph.read_ped(DATA / "nat2_trio.ped")
        ev["trio"], notes = ph.trio_orientations(vcf, CHILD, *family[CHILD])
        assert not notes
    if reads:
        ev["reads"] = ph.remap(read_vcf(DATA / "phasing" / "nat2_trio.whatshap.vcf"), vcf, CHILD)
    if panel:
        ev["vcf"] = ph.remap(panel, vcf, CHILD)
    return ev


def test_unphased_input_has_no_evidence(nat2):
    assert ph.orientations_from_vcf(nat2, CHILD) == {}
    assert audit_sample(nat2, CHILD).n_scenarios == 32  # 6 hétérozygotes non phasés


def test_trio_and_reads_agree_trans(nat2):
    fused = ph.fuse(nat2, CHILD, _evidence(nat2))
    assert fused.sources == ["trio", "reads"]
    p = _pair(fused, nat2, I114T, R268K)
    assert p.by_source == {"trio": "trans", "reads": "trans"} and p.final == "trans" and p.status == "accord"
    assert _pair(fused, nat2, I114T, L161).final == "cis"
    assert not fused.disagreements and not fused.notes


def test_fused_vcf_gives_single_scenario(nat2):
    fused = ph.fuse(nat2, CHILD, _evidence(nat2))
    audit = audit_sample(ph.apply(nat2, {CHILD: fused}), CHILD)
    assert audit.n_scenarios == 1
    assert audit.units == ["bloc 1 (trio + reads (WhatsHap))"]


def test_wrong_input_phase_is_flagged(nat2):
    panel = read_vcf(DATA / "phasing" / "nat2_trio.panel.vcf")
    i = _idx(panel, R268K)
    gt = panel.records[i].genotypes[CHILD]
    panel.records[i].genotypes[CHILD] = Genotype(gt.alleles[::-1], True, None)  # erreur de switch simulée
    fused = ph.fuse(nat2, CHILD, _evidence(nat2, reads=False, panel=panel))
    p = _pair(fused, nat2, I114T, R268K)
    assert p.by_source == {"trio": "trans", "vcf": "cis"} and p.status == "désaccord" and p.final == "trans"
    assert {(d.label_a, d.label_b) for d in fused.disagreements} >= {(p.label_a, p.label_b)}


def test_panel_phase_matches_trio(nat2):
    panel = read_vcf(DATA / "phasing" / "nat2_trio.panel.vcf")
    fused = ph.fuse(nat2, CHILD, _evidence(nat2, panel=panel))
    assert not fused.disagreements


# --- cas construits ---------------------------------------------------------------------------


def _vcf(rows, samples=("C", "P", "M")):
    """rows : (pos, {échantillon: GT}) ; GT '0/1', '0|1' ou '0|1:PS'."""
    records = []
    for pos, gts in rows:
        g = {}
        for s in samples:
            gt, _, ps = gts.get(s, "0/0").partition(":")
            sep = "|" if "|" in gt else "/"
            g[s] = Genotype(tuple(int(a) for a in gt.split(sep)), sep == "|", ps or None)
        records.append(Record("1", pos, f"v{pos}", "A", ["G"], g))
    return Vcf([], list(samples), records, True)


def test_trio_links_two_read_blocks():
    vcf = _vcf([(100, {"C": "0/1", "P": "0/1", "M": "0/1"}),   # trio non informatif
                (200, {"C": "0/1", "P": "0/0", "M": "1/1"}),   # paternel = 0
                (300, {"C": "0/1", "P": "1/1", "M": "0/0"}),   # paternel = 1
                (400, {"C": "0/1", "P": "0/1", "M": "0/1"})])
    reads = {0: ph.Orientation("b1", False), 1: ph.Orientation("b1", False),
             2: ph.Orientation("b2", False), 3: ph.Orientation("b2", True)}
    trio, _ = ph.trio_orientations(vcf, "C", "P", "M")
    assert set(trio) == {1, 2}
    fused = ph.fuse(vcf, "C", {"trio": trio, "reads": reads})
    by = {(p.a, p.b): p for p in fused.pairs}
    assert by[(1, 2)].by_source == {"trio": "trans"}
    # 100 cis 200 (reads), 200 trans 300 (trio), 300 trans 400 (reads) → 100 cis 400, déduit
    assert by[(0, 3)].by_source == {} and by[(0, 3)].status == "déduite" and by[(0, 3)].final == "cis"
    assert audit_sample(ph.apply(vcf, {"C": fused}), "C").n_scenarios == 1


def test_mendelian_error_is_reported():
    vcf = _vcf([(100, {"C": "0/1", "P": "0/0", "M": "0/0"})])
    orient, notes = ph.trio_orientations(vcf, "C", "P", "M")
    assert orient == {} and "incompatibilité mendélienne" in notes[0]


def test_indirect_contradiction_is_noted():
    vcf = _vcf([(100, {"C": "0/1"}), (200, {"C": "0/1"}), (300, {"C": "0/1"})])
    trio = {0: ph.Orientation("trio", False), 1: ph.Orientation("trio", True)}    # 100 trans 200
    reads = {1: ph.Orientation("r", False), 2: ph.Orientation("r", False)}        # 200 cis 300
    given = {0: ph.Orientation("g", False), 2: ph.Orientation("g", False)}        # 100 cis 300 : contredit
    fused = ph.fuse(vcf, "C", {"trio": trio, "reads": reads, "vcf": given})
    by = {(p.a, p.b): p for p in fused.pairs}
    assert by[(0, 2)].final == "trans" and by[(0, 2)].status == "source unique"
    assert len(fused.notes) == 1 and "contredit par d'autres relations" in fused.notes[0]


def test_read_ped(tmp_path):
    ped = tmp_path / "f.ped"
    ped.write_text("#fam ind pere mere\nF C P 0 1 0\nF P 0 0 1 0\n")
    assert ph.read_ped(ped) == {"C": ("P", None), "P": (None, None)}
