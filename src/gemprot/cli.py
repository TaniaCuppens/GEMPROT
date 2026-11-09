"""Command-line interface: `gemprot indiv …`."""

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
from .phasing import SOURCE_LABELS, read_ped
from .scope import Outside

app = typer.Typer(help="GEMPROT 2: the two protein copies of a sample, from a phased VCF.")


@app.callback()
def main() -> None:
    pass


def _print_outside(name: str, outside: list[Outside]) -> None:
    carried = [site for site in outside if site.record.genotypes[name].carried]
    if not carried:
        return
    for site in carried:
        if site.splice:
            typer.echo(f"  ! {site.record.label(1)}: {site.consequence}. The splicing effect is not modeled")
    listed = ", ".join(f"{site.record.label(1)} ({site.consequence})" for site in carried)
    typer.echo(f"  excluded from the scenarios (no change to the protein): {listed}")


def _print_sample(
    name: str, res: indiv_mod.SampleResult, ref_cds: str, ann: Annotations | None, outside: list[Outside]
) -> None:
    audit = res.audit
    typer.echo(f"{name}")
    for note in audit.notes:
        typer.echo(f"  ! {note}")
    _print_outside(name, outside)
    if res.phase:
        f = res.phase
        typer.echo(f"  phase sources: {', '.join(SOURCE_LABELS[s] for s in f.sources) or 'none'}")
        for note in f.notes:
            typer.echo(f"  ! {note}")
        for p in f.disagreements:
            opinions = ", ".join(f"{SOURCE_LABELS[s]}: {r}" for s, r in p.by_source.items())
            typer.echo(f"  ! conflict {p.label_a} / {p.label_b} ({opinions}) → kept: {p.final}")
    evidence = {(e.a, e.b): e for e in res.phase.pairs} if res.phase else {}
    for p in audit.pairs:
        e = evidence.get((min(p.a.record, p.b.record), max(p.a.record, p.b.record)))
        if e and e.by_source:
            opinions = ", ".join(f"{SOURCE_LABELS[s]}: {r}" for s, r in e.by_source.items())
            typer.echo(f"  phase {p.a.label} / {p.b.label}: {p.relation} ({opinions} — {e.status})")
        else:
            typer.echo(f"  phase {p.a.label} / {p.b.label}: {p.relation} ({p.reason})")
    many = len(res.scenarios) > 1
    if many:
        verdict = "same conclusion for all phases" if res.conclusive else "the conclusion DEPENDS on the phase"
        typer.echo(f"  {len(res.scenarios)} phase scenarios — {verdict}")
    for s in res.scenarios:
        indent = "  "
        if many:
            rel = ", ".join([f"{a} {r} {b}" for a, b, r in s.relations] + list(s.scenario.assumptions))
            typer.echo(f"  scenario {s.scenario.index + 1} ({rel})")
            indent = "    "
        typer.echo(f"{indent}→ {s.verdict.summary}")
        for i, (copy, alleles) in enumerate(zip(s.verdict.copies, s.scenario.copies), 1):
            typer.echo(f"{indent}  copy {i}: {copy.label}")
            if copy.masked:
                masked = ", ".join(hgvs_c(d, ref_cds) for d in copy.masked)
                typer.echo(f"{indent}    no effect of their own (downstream of the truncation, in cis): {masked}")
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
                typer.echo(f"{indent}   {flag} {a.label} ({', '.join(bits)}): {cv.status}{where}")


@app.command()
def indiv(
    vcf: Path = typer.Argument(..., exists=True, help="VCF, phased or unphased"),
    transcript: str = typer.Option(..., "--transcript", "-t", help="Transcript, for example NM_004004.6"),
    region: Optional[str] = typer.Option(None, "--region", "-r", help="Region to extract (indexed VCF)"),
    sample: list[str] = typer.Option([], "--sample", "-s", help="Sample to keep. You can repeat this option"),
    bam: list[Path] = typer.Option(
        [], "--bam", exists=True, help="Indexed BAM/CRAM for read phasing with WhatsHap. You can repeat this option"
    ),
    ped: Optional[Path] = typer.Option(
        None, "--ped", exists=True, help="Pedigree .ped file for trio phasing. The parents must be in the VCF"
    ),
    html: Optional[Path] = typer.Option(None, "--html", help="Write the standalone HTML report to this path"),
    annotate: bool = typer.Option(True, help="Add ClinVar, AlphaMissense, REVEL, CADD and protein domains"),
    protein_only: bool = typer.Option(
        True, "--protein-only/--all-variants",
        help="Use only the variants that change the protein in the scenarios. "
        "The report lists the other variants separately",
    ),
    offline: bool = typer.Option(False, help="Send no request to UniProt/InterPro. Use the cache only"),
    config: Path = typer.Option(config_mod.DEFAULT_CONFIG, help="Configuration file"),
) -> None:
    """For each sample: phase of the variants, protein copies with a changed sequence, variants masked in cis."""
    cfg = config_mod.load(config)
    with tempfile.TemporaryDirectory(prefix="gemprot_", dir=Path.cwd()) as tmp:
        tmpdir = Path(tmp)
        vcf_in = vcf.resolve()
        if region:
            keep = list(sample)
            if sample and ped:  # keep the parents for trio phasing
                family = read_ped(ped)
                keep += [p for s in sample for p in family.get(s, ()) if p and p not in keep]
            vcf_in = extract_region(vcf_in, region, tmpdir / "region.vcf.gz", cfg, keep)
        try:
            result = indiv_mod.run(vcf_in, transcript, cfg, tmpdir, sample, annotate, not offline,
                                   [b.resolve() for b in bam], ped.resolve() if ped else None, protein_only)
        except (KeyError, ValueError) as e:
            raise typer.BadParameter(str(e.args[0]))
    if result.annotations:
        for note in result.annotations.notes:
            typer.echo(f"! {note}")
    for name, res in result.samples.items():
        _print_sample(name, res, result.ref_cds, result.annotations, result.outside)
    if html:
        from .report import render_indiv

        render_indiv(result, html, source=vcf.name)
        typer.echo(f"report: {html}")
