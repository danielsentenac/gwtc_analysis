"""Slurm job chains: one sbatch script per step, submitted one after the other with afterok dependencies.

Used by the modes whose heavy steps run on a cluster (`galaxy_catalog`, `hubble_constant` with
``executor="slurm"``). The steps run a self-contained runner script with the Python interpreter of the analysis
environment; the runner is copied next to the scripts, so the jobs do not need gwtc_analysis installed.
On CC-IN2P3, every job must give its time limit, CPU count and memory, and jobs using /sps declare
``--licenses=sps``.
"""
from __future__ import annotations

import shlex
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

ARRAY_INDEX = "$SLURM_ARRAY_TASK_ID"


@dataclass
class Step:
    """One sbatch script: `args` of the runner (ARRAY_INDEX and other $-tokens are left unquoted), as an array job
    of `array` tasks when given, with extra #SBATCH `options` (after the common ones, so they win) and shell lines
    `pre` run before the runner."""
    name: str
    args: list[str]
    array: Optional[int] = None
    options: list[str] = field(default_factory=list)
    pre: list[str] = field(default_factory=list)


def _sbatch_lines(opts: Iterable[str]) -> str:
    return "".join(f"#SBATCH {o if o.startswith('--') else '--' + o}\n" for o in opts)


def _arg(a: str) -> str:
    return a if a.startswith("$") or a.startswith("${") else shlex.quote(a)


def write_chain(directory: Path, python: str, runner: Path, steps: list[Step], *, logs: Path,
                prefix: str, options: Iterable[str] = (), env_setup: str = "",
                env: Optional[dict] = None) -> Path:
    """Write <directory>/NN_<step>.sh for each step and <directory>/submit.sh, which submits them in order, each
    after the previous one succeeded. Returns submit.sh."""
    directory.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    local_runner = directory / runner.name
    shutil.copy2(runner, local_runner)
    common = list(options)
    names = []
    for i, st in enumerate(steps):
        log = logs / (f"{st.name}_%A_%a.log" if st.array else f"{st.name}_%j.log")
        body = f"#!/bin/bash\n#SBATCH --job-name={prefix}_{st.name}\n#SBATCH --output={log}\n"
        if st.array:
            body += f"#SBATCH --array=0-{st.array - 1}\n"
        body += _sbatch_lines(common) + _sbatch_lines(st.options)
        body += "set -euo pipefail\n" + (env_setup.rstrip() + "\n" if env_setup else "")
        body += "export HDF5_USE_FILE_LOCKING=FALSE\n"
        ld = (env or {}).get("LD_LIBRARY_PATH", "")
        if ld:
            body += f"export LD_LIBRARY_PATH={shlex.quote(ld)}${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}\n"
        body += "".join(line.rstrip() + "\n" for line in st.pre)
        body += " ".join([shlex.quote(python), shlex.quote(str(local_runner))] + [_arg(a) for a in st.args]) + "\n"
        f = directory / f"{i:02d}_{st.name}.sh"
        f.write_text(body)
        f.chmod(0o755)
        names.append(f.name)
    sub = ["#!/bin/bash", "# submits the steps, each after the previous one succeeded", "set -euo pipefail",
           f"cd {shlex.quote(str(directory))}", 'dep=""']
    for n in names:
        sub.append(f'jid=$(sbatch --parsable ${{dep:+--dependency=afterok:$dep}} {n}); echo "{n}: job $jid"; dep=$jid')
    s = directory / "submit.sh"
    s.write_text("\n".join(sub) + "\n")
    s.chmod(0o755)
    return s


def submit(script: Path) -> str:
    import subprocess

    r = subprocess.run(["bash", str(script)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"sbatch failed: {r.stderr.strip()}")
    return r.stdout
