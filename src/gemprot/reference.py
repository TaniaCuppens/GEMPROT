"""Lecture d'une FASTA indexée (.fai), sans dépendance ; accepte « 13 » ou « chr13 »."""

from __future__ import annotations

from pathlib import Path


class Fasta:
    def __init__(self, path: Path | str):
        self.path = Path(path)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with open(f"{self.path}.fai") as fh:
            for line in fh:
                name, length, offset, bases, width = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(bases), int(width))

    def _name(self, chrom: str) -> str:
        for cand in (chrom, f"chr{chrom}", chrom.removeprefix("chr")):
            if cand in self.index:
                return cand
        raise KeyError(f"chromosome {chrom} absent de {self.path}")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Séquence 1-based inclusive [start, end]."""
        length, offset, bases, width = self.index[self._name(chrom)]
        start, end = max(start, 1), min(end, length)
        first = offset + (start - 1) // bases * width + (start - 1) % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases
        with open(self.path, "rb") as fh:
            fh.seek(first)
            raw = fh.read(last - first + 1).decode()
        return raw.replace("\n", "").replace("\r", "").upper()
