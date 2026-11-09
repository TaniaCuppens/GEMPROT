"""Haplosaurus (`haplo`): run it and read its JSON output.

haplo applies all the variants of a haplotype to the coding sequence (CDS) of a transcript.
GEMPROT keeps its haplotype CDS (sequence, differences, carrier samples). Then `classify`
translates and interprets them.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .tools import run_tool


@dataclass(frozen=True)
class CdsDiff:
    diff: str  # haplo notation, in CDS coordinates: "224G>A", "30delG"
    variant_id: str | None


@dataclass
class CdsHaplotype:
    name: str
    seq: str  # haplotype CDS, without dashes
    diffs: list[CdsDiff]
    samples: dict[str, int]  # sample → number of copies (1 or 2)


@dataclass
class TranscriptHaplotypes:
    transcript_id: str
    ref_cds: str
    haplotypes: list[CdsHaplotype] = field(default_factory=list)

    def sample_haplotypes(self) -> dict[str, list[CdsHaplotype]]:
        """For each sample, its two copies (or one copy if haploid)."""
        out: dict[str, list[CdsHaplotype]] = {}
        for hap in self.haplotypes:
            for sample, n in hap.samples.items():
                out.setdefault(sample, []).extend([hap] * n)
        return out


def _ref_cds(cds_haplotypes: list[dict]) -> str:
    for h in cds_haplotypes:
        if not h.get("diffs"):
            return h["seq"]
    # no carrier of the reference: the first aligned sequence is the reference
    return cds_haplotypes[0]["aligned_sequences"][0].replace("-", "")


def parse_haplo_json(path: Path | str) -> list[TranscriptHaplotypes]:
    result = []
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            cds = rec.get("cds_haplotypes") or []
            if not cds:
                continue
            tx = TranscriptHaplotypes(rec["transcript_id"], _ref_cds(cds))
            for h in cds:
                tx.haplotypes.append(
                    CdsHaplotype(
                        name=h["name"],
                        seq=h["seq"],
                        diffs=[CdsDiff(d["diff"], d.get("variation_feature")) for d in h.get("diffs", [])],
                        samples=dict(h.get("samples", {})),
                    )
                )
            result.append(tx)
    return result


def run_haplo(vcf: Path, out_json: Path, cfg: Config) -> Path:
    vep = cfg.vep
    cmd = [
        "haplo", "--cache", "--offline",
        "--dir_cache", vep.cache_dir,
        "--cache_version", str(vep.cache_version),
        "--assembly", cfg.assembly,
        "-i", vcf, "-o", out_json, "--json", "--force_overwrite",
    ]
    if vep.cache_type == "refseq":
        cmd.insert(3, "--refseq")
    run_tool(cmd, vep.container, (*cfg.binds, vep.cache_dir, vcf.parent, out_json.parent))
    return out_json
