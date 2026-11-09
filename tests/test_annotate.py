"""Annotations from frozen sources (tests/make_annotation_fixtures.py): no container, no network."""

import json
import re
from pathlib import Path

import pytest

from gemprot import annotate as an
from gemprot import indiv
from gemprot.classify import translate
from gemprot.haplotypes import parse_haplo_json
from gemprot.report import render_indiv

DATA = Path(__file__).parent / "data"
ANNOT = DATA / "annot"


def _load(name, transcript):
    vcf, audits, plan = indiv.prepare(DATA / f"{name}.grch38.vcf")
    tx = {t.transcript_id: t for t in parse_haplo_json(DATA / "haplo" / f"{name}.scenarios.haplo.json")}[transcript]
    vep = an.parse_vep(ANNOT / f"{name}.vep.json", transcript)
    clinvar = an.parse_clinvar((ANNOT / f"{name}.clinvar.vcf").read_text().splitlines())
    lines = (ANNOT / f"{name}.dbnsfp.tsv").read_text().splitlines()
    scores = an.parse_dbnsfp(lines[0], lines[1:], vep.mane_ensembl)
    protein = an.protein_info_from_json(json.loads((ANNOT / f"{name}.protein.json").read_text()))
    length = len(translate(tx.ref_cds).rstrip("*"))
    ann = an.combine(vcf, vep, clinvar, scores, protein, length)
    return indiv.IndivResult(transcript, tx.ref_cds, indiv.assemble(tx, audits, plan), ann, vcf)


def _statuses(result, sample):
    (sc,) = result.samples[sample].scenarios
    ann = result.annotations
    return [
        [(ann.variants[cv.key].label, cv.status, cv.effect_on_copy) for cv in an.copy_variants(c, al, ann)]
        for c, al in zip(sc.verdict.copies, sc.scenario.copies)
    ]


@pytest.fixture(scope="module")
def gjb2():
    return _load("gjb2_r75q_35delG", "NM_004004.6")


@pytest.fixture(scope="module")
def synth():
    return _load("prnp_synthetic", "NM_000311.5")


def test_sources(gjb2):
    ann = gjb2.annotations
    r75q = next(v for v in ann.variants.values() if v.label == "rs28931593")
    assert r75q.isolated.hgvsp == "p.Arg75Gln" and r75q.isolated.hgvsc == "c.224G>A"
    assert r75q.clinvar.pathogenic and r75q.clinvar.stars == 2
    assert r75q.scores.alphamissense_class == "likely pathogenic"
    dfn = next(v for v in ann.variants.values() if v.label == "rs80338939")
    assert dfn.isolated.hgvsp == "p.Gly12Valfs*2" and dfn.scores is None  # dbNSFP: SNVs only
    assert ann.protein.uniprot == "P29033" and [d.accession for d in ann.protein.domains] == ["PF00029"]
    assert not ann.notes


def test_dbnsfp_uses_mane_transcript(synth):
    """D178N: one AlphaMissense score per transcript. GEMPROT keeps the MANE score."""
    d178n = next(v for v in synth.annotations.variants.values() if v.label == "c532")
    assert d178n.scores.alphamissense == pytest.approx(0.7837)


def test_masked_pathogenic_variant_is_flagged(gjb2):
    copy1, copy2 = _statuses(gjb2, "CIS")
    assert ("rs28931593", "masked", "p.Gly12Valfs*2") in copy1
    assert ("rs80338939", "applied", None) in copy1
    assert copy2 == []


def test_trans_variants_apply(gjb2):
    assert [s for copy in _statuses(gjb2, "TRANS") for _, s, _ in copy] == ["applied", "applied"]


@pytest.mark.parametrize(
    "sample, expected",
    [
        ("MEME_CODON_CIS", [("c532", "combined", "p.Asp178Lys"), ("c534", "combined", "p.Asp178Lys")]),
        ("FS_COMPENSE", [("fs500", "combined", "p.Asp167_Glu168delinsValThr"),
                         ("ins503", "combined", "p.Asp167_Glu168delinsValThr")]),
        ("FS_STOP", [("fs400", "applied", None), ("ins409", "masked", "p.Met134*")]),
        ("STOP_PERDU", [("stoplost", "applied", None)]),  # VEP: ext*16, haplo: ext*?
        ("START_PERDU", [("start", "applied", None)]),
    ],
)
def test_copy_status(synth, sample, expected):
    copy1, _ = _statuses(synth, sample)
    assert copy1 == expected


def test_same_codon_in_trans_each_applies(synth):
    copy1, copy2 = _statuses(synth, "MEME_CODON_TRANS")
    assert copy1 == [("c532", "applied", None)] and copy2 == [("c534", "applied", None)]


def test_report_with_annotations(gjb2, tmp_path):
    page = render_indiv(gjb2, tmp_path / "r.html", source="x.vcf").read_text()
    assert not re.search(r"""src=["']?(https?:)?//""", page) and "<link" not in page
    assert "UniProt P29033" in page and "Connexin" in page
    assert "p.Gly12Valfs*2 masks" in page  # alert of the CIS case
    assert page.count("<li>") == 1  # one alert only (CIS)
    assert "VCV 17027" in page and "0.87" in page
