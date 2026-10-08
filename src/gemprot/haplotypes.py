"""Haplosaurus (`haplo`) : lancement et lecture de la sortie JSON.

haplo applique tous les variants d'un haplotype à la séquence codante (CDS) d'un transcrit.
On garde ses CDS d'haplotype (séquence + différences + échantillons porteurs) ; la traduction
et l'interprétation sont faites ensuite par `classify`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config
from .tools import run_tool


@dataclass(frozen=True)
class CdsDiff:
    diff: str  # notation haplo, en coordonnées CDS : "224G>A", "30delG"
    variant_id: str | None


@dataclass
class CdsHaplotype:
    name: str
    seq: str  # CDS de l'haplotype, sans tirets
    diffs: list[CdsDiff]
    samples: dict[str, int]  # échantillon → nombre de copies (1 ou 2)


@dataclass
class TranscriptHaplotypes:
    transcript_id: str
    ref_cds: str
    haplotypes: list[CdsHaplotype] = field(default_factory=list)

    def sample_haplotypes(self) -> dict[str, list[CdsHaplotype]]:
        """Pour chaque échantillon, ses deux copies (ou une seule si haploïde)."""
        out: dict[str, list[CdsHaplotype]] = {}
        for hap in self.haplotypes:
            for sample, n in hap.samples.items():
                out.setdefault(sample, []).extend([hap] * n)
        return out


def _ref_cds(cds_haplotypes: list[dict]) -> str:
    for h in cds_haplotypes:
        if not h.get("diffs"):
            return h["seq"]
    # aucun porteur de la référence : la 1re séquence alignée est la référence
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
