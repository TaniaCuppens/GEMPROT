"""Cas du papier (PRNP, HFE) et cas synthétiques de classify, sur sorties haplo figées.

VCF : tests/make_test_vcfs.py ; sorties haplo : tests/data/haplo/<nom>.scenarios.haplo.json.
"""

from pathlib import Path

import pytest

from gemprot import indiv
from gemprot.haplotypes import parse_haplo_json

DATA = Path(__file__).parent / "data"


def _results(name, transcript):
    vcf, audits, plan = indiv.prepare(DATA / f"{name}.grch38.vcf")
    txs = {t.transcript_id: t for t in parse_haplo_json(DATA / "haplo" / f"{name}.scenarios.haplo.json")}
    return indiv.assemble(txs[transcript], audits, plan)


def _copies(results, sample):
    (scenario,) = results[sample].scenarios
    return scenario.verdict.copies


# --- PRNP (NM_000311.5, brin +) : Monari 1994 --------------------------------------------------


@pytest.fixture(scope="module")
def prnp():
    return _results("prnp_d178n_m129v", "NM_000311.5")


def test_prnp_ffi_d178n_with_129m(prnp):
    """Insomnie fatale familiale : 178N porté par la copie 129M, 129V sur l'autre copie."""
    assert [c.label for c in _copies(prnp, "IFF")] == ["p.Asp178Asn", "p.Met129Val"]


def test_prnp_cjd_d178n_with_129v(prnp):
    """Creutzfeldt-Jakob : 178N et 129V sur la même copie."""
    assert [c.label for c in _copies(prnp, "MCJ")] == ["p.Met129Val, p.Asp178Asn", "référence"]


def test_prnp_same_variants_different_proteins(prnp):
    """Mêmes génotypes variant par variant (hétérozygotes) mais protéines différentes."""
    iff = {c.protein for c in _copies(prnp, "IFF")}
    mcj = {c.protein for c in _copies(prnp, "MCJ")}
    assert iff.isdisjoint(mcj)


# --- HFE (NM_000410.4, brin +) : Uguen 2017 ----------------------------------------------------


@pytest.fixture(scope="module")
def hfe():
    return _results("hfe_c282y_h63d", "NM_000410.4")


@pytest.mark.parametrize(
    "sample, labels, altered",
    [
        ("C282Y_HOMO", ["p.Cys282Tyr", "p.Cys282Tyr"], 2),
        ("C282Y_H63D_TRANS", ["p.Cys282Tyr", "p.His63Asp"], 2),
        ("C282Y_H63D_CIS", ["p.His63Asp, p.Cys282Tyr", "référence"], 1),
    ],
)
def test_hfe(hfe, sample, labels, altered):
    copies = _copies(hfe, sample)
    assert [c.label for c in copies] == labels
    assert sum(c.altered for c in copies) == altered


# --- cas synthétiques sur PRNP -----------------------------------------------------------------


@pytest.fixture(scope="module")
def synth():
    return _results("prnp_synthetic", "NM_000311.5")


@pytest.mark.parametrize(
    "sample, label, truncated_at",
    [
        ("START_PERDU", "p.Met1?", 0),
        ("STOP_GAGNE", "p.Trp7*", 6),
        ("DEL_EN_PHASE", "p.Met129del", None),
        ("STOP_PERDU", "p.Ter254Argext*?", None),  # CDS haplo sans 3'UTR : longueur inconnue
        ("MEME_CODON_CIS", "p.Asp178Lys", None),  # GAC → AAA : ni Asn ni Glu
    ],
)
def test_single_altered_copy(synth, sample, label, truncated_at):
    altered, ref = _copies(synth, sample)
    assert altered.label == label and altered.truncated_at == truncated_at
    assert ref.label == "référence"


def test_same_codon_in_trans(synth):
    assert [c.label for c in _copies(synth, "MEME_CODON_TRANS")] == ["p.Asp178Asn", "p.Asp178Glu"]


def test_multiallelic(synth):
    assert [c.label for c in _copies(synth, "MULTI")] == ["p.Met129Val", "p.Met129Leu"]


def test_frameshift_creating_stop_masks_downstream_insertion(synth):
    altered, _ = _copies(synth, "FS_STOP")
    assert altered.label == "p.Met134*"
    assert [d.diff for d in altered.masked] == ["410insA"]
    assert not altered.compensated_frameshift  # le stop arrive avant la 2e indel


def test_compensated_frameshift(synth):
    altered, _ = _copies(synth, "FS_COMPENSE")
    assert altered.compensated_frameshift
    assert altered.label == "p.Asp167_Glu168delinsValThr"
    assert altered.truncated_at is None
