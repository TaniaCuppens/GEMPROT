"""Audit de la phase (GT + PS) et scénarios cis/trans, sur GJB2 p.Arg75Gln + c.35delG.

Les sorties haplo des VCF de scénarios sont figées dans tests/data/haplo/ ; pour les régénérer,
lancer haplo sur tests/data/haplo/*.scenarios.vcf (voir test_container_regenerates_fixture).
"""

import os
from pathlib import Path

import pytest

from gemprot import indiv, phase
from gemprot.haplotypes import parse_haplo_json

DATA = Path(__file__).parent / "data"
MANE = "NM_004004.6"
R75Q, DEL35 = "rs28931593", "rs80338939"


def _results(name):
    vcf, audits, plan = indiv.prepare(DATA / f"{name}.grch38.vcf")
    txs = {t.transcript_id: t for t in parse_haplo_json(DATA / "haplo" / f"{name}.scenarios.haplo.json")}
    return indiv.assemble(txs[MANE], audits, plan)


@pytest.fixture(scope="module")
def cases():
    return _results("gjb2_phase_cases")


def _relation(res):
    (pair,) = res.audit.pairs
    return pair.relation


def _labels(scenario_result):
    return sorted(c.label for c in scenario_result.verdict.copies)


def test_same_ps_cis(cases):
    r = cases["CIS_MEME_PS"]
    assert _relation(r) == "cis" and len(r.scenarios) == 1
    assert _labels(r.scenarios[0]) == ["p.Gly12Valfs*2", "référence"]
    assert r.scenarios[0].verdict.altered_copies == 1


def test_same_ps_trans(cases):
    r = cases["TRANS_MEME_PS"]
    assert _relation(r) == "trans"
    assert _labels(r.scenarios[0]) == ["p.Arg75Gln", "p.Gly12Valfs*2"]


def test_copy_order_follows_vcf(cases):
    """Copie 1 = gauche du '|' : TRANS_MEME_PS est 0|1 pour R75Q et 1|0 pour c.35delG."""
    copies = cases["TRANS_MEME_PS"].scenarios[0].verdict.copies
    assert [c.label for c in copies] == ["p.Gly12Valfs*2", "p.Arg75Gln"]


@pytest.mark.parametrize("sample", ["PS_DIFFERENTS", "NON_PHASE"])
def test_unknown_phase_gives_both_scenarios(cases, sample):
    r = cases[sample]
    assert _relation(r) == "inconnue"
    by_relation = {s.relations[0][2]: s for s in r.scenarios}
    assert set(by_relation) == {"cis", "trans"}
    assert by_relation["cis"].verdict.altered_copies == 1
    assert by_relation["trans"].verdict.altered_copies == 2
    assert not r.conclusive  # le diagnostic dépend de la phase : il faut le dire


def test_homozygous(cases):
    r = cases["HOMOZYGOTE"]
    assert _relation(r) == "homozygote" and len(r.scenarios) == 1
    assert _labels(r.scenarios[0]) == ["p.Arg75Gln", "p.Gly12Valfs*2"]


def test_partially_missing_gives_both_possibilities(cases):
    """R75Q 1|0 et c.35delG 0|. : la copie 2 peut porter ou non c.35delG."""
    r = cases["MANQUANT"]
    assert r.audit.pairs == []  # c.35delG n'est porté par aucune copie connue
    assert any(DEL35 in n and "partiellement manquant" in n for n in r.audit.notes)
    by_assumption = {s.scenario.assumptions: s for s in r.scenarios}
    absent = by_assumption[(f"{DEL35} absent sur la copie 2",)]
    present = by_assumption[(f"{DEL35} présent sur la copie 2",)]
    assert _labels(absent) == ["p.Arg75Gln", "référence"]
    assert _labels(present) == ["p.Arg75Gln", "p.Gly12Valfs*2"]
    assert not r.conclusive


def test_haploid(cases):
    r = cases["HAPLOIDE"]
    assert _relation(r) == "cis" and r.audit.ploidy == 1
    (copy,) = r.scenarios[0].verdict.copies
    assert copy.label == "p.Gly12Valfs*2" and [d.diff for d in copy.masked] == ["224G>A"]
    assert r.scenarios[0].verdict.summary == "la copie altérée"


def test_vcf_without_ps_is_globally_phased():
    res = _results("gjb2_r75q_35delG")
    assert _relation(res["CIS"]) == "cis"
    assert _relation(res["TRANS"]) == "trans"
    assert "supposée globale" in res["CIS"].audit.units[0]


def test_too_many_scenarios():
    calls = phase.SampleAudit("X", 2, [], [f"u{i}" for i in range(6)], [])
    with pytest.raises(ValueError, match="il faut phaser"):
        phase.scenarios(calls)


def test_pseudo_name_roundtrip():
    assert phase.parse_pseudo_name(phase.pseudo_name("NA12878::x", 3, 1)) == ("NA12878::x", 3, 1)


@pytest.mark.container
@pytest.mark.skipif(not os.environ.get("GEMPROT_CONTAINER_TESTS"), reason="GEMPROT_CONTAINER_TESTS non défini")
def test_container_regenerates_fixture(tmp_path):
    """Chaîne complète avec haplo : même résultat que la sortie figée."""
    from gemprot import config

    live = indiv.run(DATA / "gjb2_phase_cases.grch38.vcf", MANE, config.load(), tmp_path)
    frozen = _results("gjb2_phase_cases")
    for sample, res in frozen.items():
        assert [_labels(s) for s in live.samples[sample].scenarios] == [_labels(s) for s in res.scenarios]
