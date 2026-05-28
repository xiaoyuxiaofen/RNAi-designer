from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from .accessibility import rnaplfold_available
from .bowtie import resolve_bowtie_build_executable, resolve_bowtie_executable
from .bowtie2 import resolve_bowtie2_build_executable, resolve_bowtie2_executable


@dataclass(frozen=True)
class DependencyStatus:
    name: str
    available: bool
    detail: str


def check_dependencies() -> list[DependencyStatus]:
    return [
        _python_status(),
        _viennarna_python_status(),
        _tool_status("RNAplfold executable", rnaplfold_available("RNAplfold"), shutil.which("RNAplfold") or ""),
        _tool_status("Bowtie 1 aligner", bool(resolve_bowtie_executable("bowtie")), resolve_bowtie_executable("bowtie") or ""),
        _tool_status(
            "Bowtie 1 indexer",
            bool(resolve_bowtie_build_executable("bowtie-build")),
            resolve_bowtie_build_executable("bowtie-build") or "",
        ),
        _tool_status("Bowtie2 aligner", bool(resolve_bowtie2_executable("bowtie2")), resolve_bowtie2_executable("bowtie2") or ""),
        _tool_status(
            "Bowtie2 indexer",
            bool(resolve_bowtie2_build_executable("bowtie2-build")),
            resolve_bowtie2_build_executable("bowtie2-build") or "",
        ),
    ]


def format_dependency_statuses(statuses: list[DependencyStatus]) -> str:
    lines = ["Dependency check:"]
    for status in statuses:
        marker = "OK" if status.available else "MISSING"
        detail = f" - {status.detail}" if status.detail else ""
        lines.append(f"[{marker}] {status.name}{detail}")
    return "\n".join(lines)


def _python_status() -> DependencyStatus:
    return DependencyStatus(
        name="Python",
        available=True,
        detail=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro} at {sys.executable}",
    )


def _viennarna_python_status() -> DependencyStatus:
    _add_local_python_packages()
    try:
        import RNA  # type: ignore[import-not-found]
    except ImportError:
        return DependencyStatus("ViennaRNA Python bindings", False, "")
    return DependencyStatus("ViennaRNA Python bindings", True, f"version {getattr(RNA, '__version__', 'unknown')}")


def _tool_status(name: str, available: bool, detail: str) -> DependencyStatus:
    return DependencyStatus(name=name, available=available, detail=detail)


def _add_local_python_packages() -> None:
    project_root = Path(__file__).resolve().parents[2]
    local_packages = project_root / "tools" / "python-packages"
    if local_packages.exists():
        local_path = str(local_packages)
        if local_path not in sys.path:
            sys.path.insert(0, local_path)
