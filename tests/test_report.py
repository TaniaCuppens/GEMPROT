"""Rapport HTML : autonome (aucune ressource externe) et fidèle aux résultats."""

import re
from pathlib import Path

from gemprot import indiv
from gemprot.haplotypes import parse_haplo_json
from gemprot.report import render_indiv

DATA = Path(__file__).parent / "data"


def _render(tmp_path):
    vcf, audits, plan = indiv.prepare(DATA / "gjb2_phase_cases.grch38.vcf")
    tx = {t.transcript_id: t for t in parse_haplo_json(DATA / "haplo" / "gjb2_phase_cases.scenarios.haplo.json")}["NM_004004.6"]
    result = indiv.IndivResult("NM_004004.6", tx.ref_cds, indiv.assemble(tx, audits, plan))
    return render_indiv(result, tmp_path / "r.html", source="gjb2_phase_cases.grch38.vcf").read_text()


def test_report_is_self_contained(tmp_path):
    page = _render(tmp_path)
    assert not re.search(r"""src=["']?(https?:)?//""", page) and "<link" not in page
    assert "@import" not in page and "url(http" not in page


def test_report_content(tmp_path):
    page = _render(tmp_path)
    assert "GEMPROT · NM_004004.6" in page and "226 acides aminés" in page
    assert page.count('class="card" id="s-') == 7  # un bloc par individu
    assert "Conclusion dépendante de la phase" in page  # PS_DIFFERENTS, NON_PHASE, MANQUANT
    assert "Gly12Valfs*2" in page and "c.224G&gt;A" in page
    assert "sans effet propre" in page  # info-bulle du variant masqué
    assert page.count("Scénario 2") == 3
