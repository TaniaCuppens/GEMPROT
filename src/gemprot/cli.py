"""Interface en ligne de commande : `gemprot indiv …`."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import typer

from . import config as config_mod
from . import indiv as indiv_mod
from .annotate import Annotations, copy_variants
from .classify import hgvs_c
from .extract import extract_region

app = typer.Typer(help="GEMPROT 2 : les deux protéines d'un individu à partir d'un VCF phasé.")


@app.callback()
def main() -> None:
    pass


def _print_sample(name: str, res: indiv_mod.SampleResult, ref_cds: str, ann: Annotations | None) -> None:
    audit = res.audit
    typer.echo(f"{name}")
    for note in audit.notes:
        typer.echo(f"  ! {note}")
    for p in audit.pairs:
        typer.echo(f"  phase {p.a.label} / {p.b.label} : {p.relation} ({p.reason})")
    many = len(res.scenarios) > 1
    if many:
        verdict = "même conclusion quelle que soit la phase" if res.conclusive else "la conclusion DÉPEND de la phase"
        typer.echo(f"  {len(res.scenarios)} scénarios de phase — {verdict}")
    for s in res.scenarios:
        indent = "  "
        if many:
            rel = ", ".join([f"{a} {r} {b}" for a, b, r in s.relations] + list(s.scenario.assumptions))
            typer.echo(f"  scénario {s.scenario.index + 1} ({rel})")
            indent = "    "
        typer.echo(f"{indent}→ {s.verdict.summary}")
        for i, (copy, alleles) in enumerate(zip(s.verdict.copies, s.scenario.copies), 1):
            typer.echo(f"{indent}  copie {i} : {copy.label}")
            if copy.masked:
                masked = ", ".join(hgvs_c(d, ref_cds) for d in copy.masked)
                typer.echo(f"{indent}    sans effet propre (en aval de la troncature, en cis) : {masked}")
            if ann is None:
                continue
            for cv in copy_variants(copy, alleles, ann):
                a = ann.variants[cv.key]
                bits = [a.isolated.hgvsp or a.isolated.consequence if a.isolated else "?"]
                if a.clinvar:
                    bits.append(f"ClinVar {a.clinvar.significance}")
                if a.scores and a.scores.alphamissense is not None:
                    bits.append(f"AlphaMissense {a.scores.alphamissense:.2f}")
                flag = "!" if cv.misleading else " "
                where = f" → {cv.effect_on_copy}" if cv.effect_on_copy else ""
                typer.echo(f"{indent}   {flag} {a.label} ({', '.join(bits)}) : {cv.status}{where}")


@app.command()
def indiv(
    vcf: Path = typer.Argument(..., exists=True, help="VCF (phasé ou non)"),
    transcript: str = typer.Option(..., "--transcript", "-t", help="Transcrit, ex. NM_004004.6"),
    region: Optional[str] = typer.Option(None, "--region", "-r", help="Région à extraire (VCF indexé)"),
    sample: list[str] = typer.Option([], "--sample", "-s", help="Échantillon(s) à garder"),
    html: Optional[Path] = typer.Option(None, "--html", help="Écrit le rapport HTML autonome ici"),
    annotate: bool = typer.Option(True, help="ClinVar, AlphaMissense, REVEL, CADD et domaines"),
    offline: bool = typer.Option(False, help="Pas de requête UniProt/InterPro (cache seulement)"),
    config: Path = typer.Option(config_mod.DEFAULT_CONFIG, help="Fichier de configuration"),
) -> None:
    """Pour chaque individu : phase des variants, copies de la protéine altérées, variants masqués en cis."""
    cfg = config_mod.load(config)
    with tempfile.TemporaryDirectory(prefix="gemprot_", dir=Path.cwd()) as tmp:
        tmpdir = Path(tmp)
        vcf_in = vcf.resolve()
        if region:
            vcf_in = extract_region(vcf_in, region, tmpdir / "region.vcf.gz", cfg, sample)
        try:
            result = indiv_mod.run(vcf_in, transcript, cfg, tmpdir, sample, annotate, not offline)
        except (KeyError, ValueError) as e:
            raise typer.BadParameter(str(e.args[0]))
    if result.annotations:
        for note in result.annotations.notes:
            typer.echo(f"! {note}")
    for name, res in result.samples.items():
        _print_sample(name, res, result.ref_cds, result.annotations)
    if html:
        from .report import render_indiv

        render_indiv(result, html, source=vcf.name)
        typer.echo(f"rapport : {html}")
