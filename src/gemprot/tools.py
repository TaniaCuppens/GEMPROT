"""Runs the external tools (haplo, vep, bcftools…).

In the gemprot image, the tools are on the PATH and GEMPROT calls them directly.
Else (development), GEMPROT calls them with `singularity exec` on the container from the config.
The commands are always lists. GEMPROT never uses a shell.
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
        raise ToolError(f"{cmd[0]} not found, and singularity or apptainer is not installed")
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
        raise ToolError(f"{cmd[0]} failed (code {proc.returncode}):\n{proc.stderr.strip()}")
    return proc
