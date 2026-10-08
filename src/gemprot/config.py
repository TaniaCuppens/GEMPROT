"""Lecture de config/gemprot.yaml ; les chemins relatifs sont résolus depuis la racine du projet."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "gemprot.yaml"


@dataclass(frozen=True)
class VepConfig:
    container: Path | None
    cache_dir: Path
    cache_type: str
    cache_version: int


@dataclass(frozen=True)
class Config:
    assembly: str
    vep: VepConfig
    fasta: Path | None
    annotations: dict[str, Path]
    binds: tuple[Path, ...]  # dossiers à monter dans le conteneur


def _path(value: str | None, root: Path) -> Path | None:
    if not value:
        return None
    p = Path(value).expanduser()
    return p if p.is_absolute() else root / p


def load(path: Path | str = DEFAULT_CONFIG, root: Path = PROJECT_ROOT) -> Config:
    raw = yaml.safe_load(Path(path).read_text())
    v = raw["vep"]
    vep = VepConfig(
        container=_path(v.get("container"), root),
        cache_dir=_path(v["cache_dir"], root),
        cache_type=v.get("cache_type", "refseq"),
        cache_version=int(v["cache_version"]),
    )
    annotations = {k: _path(val, root) for k, val in (raw.get("annotations") or {}).items()}
    binds = [_path(b, root) for b in raw.get("binds", [])]
    if raw.get("genomic_dir"):
        binds.append(_path(raw["genomic_dir"], root))
    return Config(
        assembly=raw.get("assembly", "GRCh38"),
        vep=vep,
        fasta=_path(raw.get("fasta"), root),
        annotations=annotations,
        binds=tuple(binds),
    )
