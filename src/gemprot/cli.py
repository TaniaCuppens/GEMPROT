"""Interface en ligne de commande : `gemprot indiv …`."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Optional

import typer

from . import config as config_mod
from . import indiv as indiv_mod
from .classify import hgvs_c
from .extract import extract_region

app = typer.Typer(help="GEMPROT 2 : les deux protéines d'un individu à partir d'un VCF phasé.")


@app.callback()
def main() -> None:
    pass


def _print_sample(name: str, res: indiv_mod.SampleResult, ref_cds: str) -> None:
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
        for i, copy in enumerate(s.verdict.copies, 1):
            typer.echo(f"{indent}  copie {i} : {copy.label}")
            if copy.masked:
                masked = ", ".join(hgvs_c(d, ref_cds) for d in copy.masked)
                typer.echo(f"{indent}    sans effet propre (en aval de la troncature, en cis) : {masked}")


@app.command()
def indiv(
    vcf: Path = typer.Argument(..., exists=True, help="VCF (phasé ou non)"),
    transcript: str = typer.Option(..., "--transcript", "-t", help="Transcrit, ex. NM_004004.6"),
    region: Optional[str] = typer.Option(None, "--region", "-r", help="Région à extraire (VCF indexé)"),
    sample: list[str] = typer.Option([], "--sample", "-s", help="Échantillon(s) à garder"),
    html: Optional[Path] = typer.Option(None, "--html", help="Écrit le rapport HTML autonome ici"),
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
            result = indiv_mod.run(vcf_in, transcript, cfg, tmpdir, sample)
        except (KeyError, ValueError) as e:
            raise typer.BadParameter(str(e.args[0]))
    for name, res in result.samples.items():
        _print_sample(name, res, result.ref_cds)
    if html:
        from .report import render_indiv

        render_indiv(result, html, source=vcf.name)
        typer.echo(f"rapport : {html}")
