"""Extracts a region and samples from a VCF with bcftools (no shell)."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from .config import Config
from .tools import run_tool


def extract_region(
    vcf: Path,
    region: str,
    out_vcf: Path,
    cfg: Config,
    samples: Sequence[str] = (),
) -> Path:
    """Keeps the region (e.g. 13:20187470-20192974) and, if given, the requested samples.

    The input VCF must be compressed (bgzip) and indexed.
    """
    cmd = ["bcftools", "view", "-r", region, "-Oz", "-o", out_vcf]
    if samples:
        cmd += ["-s", ",".join(samples)]
    cmd.append(vcf)
    run_tool(cmd, cfg.vep.container, (*cfg.binds, vcf.parent, out_vcf.parent))
    run_tool(["bcftools", "index", "-t", out_vcf], cfg.vep.container, (*cfg.binds, out_vcf.parent))
    return out_vcf
