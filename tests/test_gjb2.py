"""Paper case: GJB2 p.Arg75Gln (dominant) + c.35delG (Iossa 2010).

In cis, c.35delG truncates the protein before codon 75. R75Q has no effect of its own, and the other
copy is normal. In trans, both copies have a changed sequence.
"""

from pathlib import Path

import pytest

from gemprot.classify import classify_transcript, hgvs_c, translate
from gemprot.haplotypes import parse_haplo_json

DATA = Path(__file__).parent / "data"
MANE = "NM_004004.6"


@pytest.fixture(scope="module")
def verdicts():
    txs = {t.transcript_id: t for t in parse_haplo_json(DATA / "haplo" / "gjb2_r75q_35delG.haplo.json")}
    return classify_transcript(txs[MANE])


@pytest.fixture(scope="module")
def ref_cds():
    txs = {t.transcript_id: t for t in parse_haplo_json(DATA / "haplo" / "gjb2_r75q_35delG.haplo.json")}
    return txs[MANE].ref_cds


def test_reference_protein(ref_cds):
    prot = translate(ref_cds)
    assert prot.startswith("MDWGTLQTILGGVNK") and prot.endswith("*")
    assert len(prot) == 227  # 226 amino acids + stop (NP_003995.2)


def test_35delG_hgvs_is_3prime_shifted(verdicts, ref_cds):
    trunc = next(c for c in verdicts["TRANS"].copies if c.truncated_at)
    delg = next(d for d in trunc.haplotype.diffs if "del" in d.diff)
    assert delg.diff == "30delG"  # haplo aligns to the left
    assert hgvs_c(delg, ref_cds) == "c.35del"


def test_cis_masks_r75q(verdicts):
    v = verdicts["CIS"]
    labels = sorted(c.label for c in v.copies)
    assert labels == ["p.Gly12Valfs*2", "reference"]
    trunc = next(c for c in v.copies if c.altered)
    assert trunc.truncated_at == 11
    assert [d.diff for d in trunc.masked] == ["224G>A"]  # R75Q masked in cis
    assert v.altered_copies == 1


def test_trans_alters_both_copies(verdicts):
    v = verdicts["TRANS"]
    assert sorted(c.label for c in v.copies) == ["p.Arg75Gln", "p.Gly12Valfs*2"]
    assert all(not c.masked for c in v.copies)
    assert v.altered_copies == 2


def test_r75q_alone(verdicts):
    v = verdicts["R75Q_SEUL"]
    assert sorted(c.label for c in v.copies) == ["p.Arg75Gln", "reference"]
    assert v.altered_copies == 1
