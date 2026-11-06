"""Paths given by jobs and by the GUI: one rule for both (swingersim.worker, swingersim.gui.server)."""

from __future__ import annotations

from pathlib import Path


def inside_root(root: Path, rel: str) -> Path:
    """Absolute path of rel under root. Refused (ValueError) when empty, when it leaves root, or
    when it contains '..'. Symlinks inside the tree are followed, not resolved: an output or
    data folder of the SWINGERsim folder may link to a disk elsewhere."""
    if rel is None or str(rel) == "":
        raise ValueError("missing path")
    if ".." in Path(str(rel)).parts:
        raise ValueError("'..' is not allowed in paths")
    p = (root / str(rel)).absolute()
    try:
        p.relative_to(root.absolute())
    except ValueError as exc:
        raise ValueError("outside the SWINGERsim folder") from exc
    return p


__all__ = ["inside_root"]
