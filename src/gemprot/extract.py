"""Extraction d'une région et d'échantillons d'un VCF avec bcftools (sans shell)."""

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
    """Garde la région (ex. 13:20187470-20192974) et, si donnés, les échantillons demandés.

    Le VCF d'entrée doit être compressé (bgzip) et indexé.
    """
    cmd = ["bcftools", "view", "-r", region, "-Oz", "-o", out_vcf]
    if samples:
        cmd += ["-s", ",".join(samples)]
    cmd.append(vcf)
    run_tool(cmd, cfg.vep.container, (*cfg.binds, vcf.parent, out_vcf.parent))
    run_tool(["bcftools", "index", "-t", out_vcf], cfg.vep.container, (*cfg.binds, out_vcf.parent))
    return out_vcf
