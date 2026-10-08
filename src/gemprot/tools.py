"""Lancement des outils externes (haplo, vep, bcftools…).

Dans l'image gemprot, les outils sont sur le PATH et sont appelés directement.
Sinon (développement), ils passent par `singularity exec` sur le conteneur de la config.
Les commandes sont toujours des listes : jamais de shell.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Sequence


class ToolError(RuntimeError):
    pass


def build_command(
    cmd: Sequence[str],
    container: Path | None = None,
    binds: Sequence[Path] = (),
) -> list[str]:
    cmd = [str(c) for c in cmd]
    if shutil.which(cmd[0]) or container is None:
        return cmd
    runner = shutil.which("singularity") or shutil.which("apptainer")
    if runner is None:
        raise ToolError(f"{cmd[0]} introuvable et ni singularity ni apptainer n'est installé")
    mounts: list[str] = []
    for b in dict.fromkeys(Path(b).resolve() for b in (*binds, Path.cwd())):
        mounts += ["-B", str(b)]
    return [runner, "exec", *mounts, str(container), *cmd]


def run_tool(
    cmd: Sequence[str],
    container: Path | None = None,
    binds: Sequence[Path] = (),
) -> subprocess.CompletedProcess[str]:
    full = build_command(cmd, container, binds)
    proc = subprocess.run(full, capture_output=True, text=True)
    if proc.returncode != 0:
        raise ToolError(f"échec de {cmd[0]} (code {proc.returncode}) :\n{proc.stderr.strip()}")
    return proc
