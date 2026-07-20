from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
from typing import Any, Mapping

PACKAGE_MODE_ALIASES = {
    "lsf": "joey_lsf",
}

PACKAGED_PROFILE_FILENAME = "hpc_profile.json"


def _clean_str(value: Any, default: str = "") -> str:
    if value is None:
        return default
    return str(value).strip()


def _clean_int(value: Any, default: int) -> int:
    try:
        return int(str(value).strip())
    except Exception:
        return int(default)


def _clean_bool(value: Any, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on", "y"}


def normalize_package_mode(
    values: Mapping[str, Any] | None,
    default_mode: str = "portable",
    lsf_enabled: bool = True,
) -> str:
    raw_default = PACKAGE_MODE_ALIASES.get((default_mode or "").strip().lower(), (default_mode or "").strip().lower())
    default = raw_default if raw_default in {"portable", "joey_lsf", "mainak_lsf", "custom_lsf"} else "portable"

    values = values or {}
    requested = PACKAGE_MODE_ALIASES.get(_clean_str(values.get("package_mode")).lower(), _clean_str(values.get("package_mode")).lower())
    include_lsf = _clean_str(values.get("include_lsf")).lower()

    if requested in {"portable", "joey_lsf", "mainak_lsf", "custom_lsf"}:
        mode = requested
    elif include_lsf in {"1", "true", "yes", "on"}:
        mode = "joey_lsf"
    elif include_lsf in {"0", "false", "no", "off"}:
        mode = "portable"
    else:
        mode = default

    if mode in {"joey_lsf", "mainak_lsf", "custom_lsf"} and not lsf_enabled:
        return "portable"
    return mode


def split_setup_commands(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item).rstrip() for item in value if str(item).strip()]
    text = str(value).replace("\r\n", "\n")
    return [line.rstrip() for line in text.split("\n") if line.strip()]


@dataclass(frozen=True)
class HPCProfile:
    profile_name: str
    scheduler: str = "lsf"
    email: str = ""
    notify_begin: bool = False
    notify_end: bool = False
    queue: str = "general"
    project: str = ""
    workers: int = 16
    mem_per_core_mb: int = 2000
    confgen_walltime: str = "48:00"
    vina_walltime: str = "96:00"
    conda_sh: str = ""
    conda_env: str = ""
    vina_executable: str = ""
    python_command: str = '$(command -v python3 || command -v python)'
    setup_commands: tuple[str, ...] = ()
    span_hosts: int = 1

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["setup_commands"] = list(self.setup_commands)
        return payload


JOEY_LSF_PROFILE = HPCProfile(
    profile_name="joey_miami",
    email="jxs794@miami.edu",
    notify_begin=True,
    notify_end=True,
    queue="general",
    project="brd",
    workers=16,
    mem_per_core_mb=2000,
    confgen_walltime="48:00",
    vina_walltime="96:00",
    conda_sh="/nethome/jxs794/miniconda3/etc/profile.d/conda.sh",
    conda_env="vina_env",
    vina_executable="$HOME/miniconda3/envs/vina_env/bin/vina",
    python_command='$(command -v python3 || command -v python)',
)

MAINAK_LSF_PROFILE = HPCProfile(
    profile_name="mainak_miami",
    email="mxb2638@miami.edu",
    notify_begin=True,
    notify_end=True,
    queue="general",
    project="brd",
    workers=16,
    mem_per_core_mb=2000,
    confgen_walltime="48:00",
    vina_walltime="96:00",
    conda_sh="/nethome/mxb2638/miniconda3/etc/profile.d/conda.sh",
    conda_env="vina_env",
    vina_executable="$HOME/miniconda3/envs/vina_env/bin/vina",
    python_command='$(command -v python3 || command -v python)',
)


def profile_from_dict(payload: Mapping[str, Any]) -> HPCProfile:
    known = {field.name for field in fields(HPCProfile)}
    setup_commands = tuple(split_setup_commands(payload.get("setup_commands")))
    cleaned = {
        "profile_name": _clean_str(payload.get("profile_name"), "custom_lsf"),
        "scheduler": _clean_str(payload.get("scheduler"), "lsf") or "lsf",
        "email": _clean_str(payload.get("email")),
        "notify_begin": _clean_bool(payload.get("notify_begin"), False),
        "notify_end": _clean_bool(payload.get("notify_end"), False),
        "queue": _clean_str(payload.get("queue"), "general") or "general",
        "project": _clean_str(payload.get("project")),
        "workers": _clean_int(payload.get("workers"), 16),
        "mem_per_core_mb": _clean_int(payload.get("mem_per_core_mb"), 2000),
        "confgen_walltime": _clean_str(payload.get("confgen_walltime"), "48:00") or "48:00",
        "vina_walltime": _clean_str(payload.get("vina_walltime"), "96:00") or "96:00",
        "conda_sh": _clean_str(payload.get("conda_sh")),
        "conda_env": _clean_str(payload.get("conda_env")),
        "vina_executable": _clean_str(payload.get("vina_executable")),
        "python_command": _clean_str(payload.get("python_command"), '$(command -v python3 || command -v python)') or '$(command -v python3 || command -v python)',
        "setup_commands": setup_commands,
        "span_hosts": _clean_int(payload.get("span_hosts"), 1),
    }
    extra = {key: value for key, value in payload.items() if key in known and key not in cleaned}
    return HPCProfile(**(cleaned | extra))


def build_custom_profile(values: Mapping[str, Any] | None) -> HPCProfile:
    values = values or {}
    return profile_from_dict(
        {
            "profile_name": _clean_str(values.get("profile_name"), "custom_lsf") or "custom_lsf",
            "email": _clean_str(values.get("lsf_email"), _clean_str(values.get("email"))),
            "notify_begin": values.get("notify_begin"),
            "notify_end": values.get("notify_end"),
            "queue": _clean_str(values.get("queue"), "general"),
            "project": _clean_str(values.get("project")),
            "workers": values.get("workers") or values.get("ncores") or 16,
            "mem_per_core_mb": values.get("mem_per_core") or values.get("mem") or 2000,
            "confgen_walltime": _clean_str(values.get("confgen_walltime"), "48:00"),
            "vina_walltime": _clean_str(values.get("vina_walltime"), _clean_str(values.get("walltime"), "96:00")),
            "conda_sh": _clean_str(values.get("conda_sh")),
            "conda_env": _clean_str(values.get("conda_env")),
            "vina_executable": _clean_str(values.get("vina_path"), _clean_str(values.get("vina_exe"))),
            "python_command": _clean_str(values.get("python_command"), '$(command -v python3 || command -v python)'),
            "setup_commands": split_setup_commands(values.get("setup_commands") or values.get("env_line")),
            "span_hosts": values.get("span_hosts") or 1,
        }
    )


def profile_for_mode(mode: str, values: Mapping[str, Any] | None = None) -> HPCProfile | None:
    if mode == "portable":
        return None
    if mode == "mainak_lsf":
        return MAINAK_LSF_PROFILE
    if mode == "custom_lsf":
        return build_custom_profile(values)
    return JOEY_LSF_PROFILE


def render_setup_block(profile: HPCProfile) -> str:
    lines: list[str] = []
    if profile.conda_sh:
        lines.append(f'source "{profile.conda_sh}"')
    if profile.conda_env:
        lines.append(f"conda activate {profile.conda_env}")
    lines.extend(profile.setup_commands)
    if not lines:
        return ""
    return "\n".join(lines) + "\n"


def render_lsf_header(
    *,
    profile: HPCProfile,
    jobname: str,
    log_prefix: str,
    walltime: str,
    workers: int | None = None,
    mem_per_core_mb: int | None = None,
) -> str:
    workers = int(workers or profile.workers)
    mem_per_core_mb = int(mem_per_core_mb or profile.mem_per_core_mb)
    lines = [
        "#!/bin/bash",
        f"#BSUB -J {jobname}",
    ]
    if profile.project:
        lines.append(f"#BSUB -P {profile.project}")
    lines.extend(
        [
            f"#BSUB -o logs/{log_prefix}_%J.out",
            f"#BSUB -e logs/{log_prefix}_%J.err",
            f"#BSUB -W {walltime}",
            f"#BSUB -q {profile.queue}",
            f"#BSUB -n {workers}",
            f'#BSUB -R "span[hosts={max(1, int(profile.span_hosts))}]"',
            f'#BSUB -R "rusage[mem={mem_per_core_mb}]"',
        ]
    )
    if profile.email:
        if profile.notify_begin:
            lines.append("#BSUB -B")
        if profile.notify_end:
            lines.append("#BSUB -N")
        if profile.notify_begin or profile.notify_end:
            lines.append(f"#BSUB -u {profile.email}")
    lines.extend(
        [
            "",
            "set -euo pipefail",
            'cd "$LS_SUBCWD"',
            'echo "PWD: $(pwd)"',
            "mkdir -p logs",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def save_packaged_profile(jobroot: Path, profile: HPCProfile) -> Path:
    path = Path(jobroot) / PACKAGED_PROFILE_FILENAME
    path.write_text(json.dumps(profile.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def load_packaged_profile(base_dir: Path | None = None) -> HPCProfile | None:
    root = Path(base_dir) if base_dir else Path.cwd()
    path = root / PACKAGED_PROFILE_FILENAME
    if not path.exists():
        return None
    return profile_from_dict(json.loads(path.read_text(encoding="utf-8")))


def packaged_profile_or_default(base_dir: Path | None = None) -> HPCProfile:
    return load_packaged_profile(base_dir) or JOEY_LSF_PROFILE


def replace_profile(profile: HPCProfile, **changes: Any) -> HPCProfile:
    if "setup_commands" in changes:
        changes["setup_commands"] = tuple(split_setup_commands(changes["setup_commands"]))
    return replace(profile, **changes)
