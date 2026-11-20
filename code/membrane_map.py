"""Pluggable membrane-plane provider for the M sub-score.

The reference plane comes from **OPM** (Orientation of Proteins in Membranes)
for experimental structures; for predicted structures the reference plane is
transferred by TM-helix superposition (Track C owns that step).  This module
defines the tiny interface MAES consumes, plus three offline implementations so
that v1.0 is runnable without network access:

* :class:`TableMembraneMap`   — load an OPM-style table (CSV/JSON) exported by Track C
* :class:`GeometricMembraneMap` — geometric fallback from TM-helix axes
* :class:`NullMembraneMap`    — no plane; the M sub-score becomes neutral (1.0)

``OPMMembraneMap.fetch`` is deliberately a stub: OPM has no public batch API,
so Track C exports the planes once and MAES loads them.

All planes use the same convention throughout MAES:
``z = (x - center) . normal`` in Angstrom, ``normal`` unit length.
"""

from __future__ import annotations

import csv
import json
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

import numpy as np
import pandas as pd

__all__ = [
    "MembranePlane",
    "MembraneMap",
    "NullMembraneMap",
    "TableMembraneMap",
    "OPMMembraneMap",
    "GeometricMembraneMap",
    "TrackBMembraneMap",
    "FrameTransform",
    "build_membrane_map",
    "load_membrane_map",
    "load_c_line_planes",
    "load_planes_auto",
    "load_frame_transforms",
    "transform_plane_opm_to_rcsb",
    "build_trackb_membrane_map",
]


def _as_vec3(value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=float).reshape(-1)
    if arr.size != 3:
        raise ValueError(f"expected 3 numbers, got {arr.size}")
    return arr


def _normalise_frame(frame: Any) -> str:
    """Normalise C-line / B-line frame labels to ``opm`` or ``rcsb``."""
    text = str(frame).strip().lower()
    if text in {"opm", "opm_oriented", "opm-oriented", "oriented"}:
        return "opm"
    if text in {"rcsb", "rcsb_original", "rcsb-original", "original", "pdb"}:
        return "rcsb"
    # tolerate decorated labels such as 'rcsb_original(mmcif)'
    if text.startswith("opm"):
        return "opm"
    if text.startswith("rcsb") or text.startswith("original"):
        return "rcsb"
    return text or "opm"


@dataclass(frozen=True)
class MembranePlane:
    """An oriented membrane plane.

    Parameters
    ----------
    normal : sequence of 3 floats
        Unit membrane normal (arbitrary global sign).
    center : sequence of 3 floats
        A point on the membrane mid-plane.
    thickness : float
        Hydrophobic-core thickness in Angstrom (OPM uses ~30 A).
    source : str
        Provenance tag, e.g. ``"opm"`` / ``"geometric"`` / ``"reference-transfer"``.
    """

    normal: Any
    center: Any
    thickness: float = 30.0
    source: str = "unknown"
    confidence: str = "high"   # high | medium | none (C-line W3 `membrane_confidence`)
    method: str = ""           # e.g. OPM_PPM / geometric / unavailable
    coord_frame: str = "opm"   # 'opm' (OPM_oriented) | 'rcsb' (RCSB_original)

    def __post_init__(self) -> None:
        normal = _as_vec3(self.normal)
        norm = float(np.linalg.norm(normal))
        if not np.isfinite(norm) or norm == 0.0:
            raise ValueError("membrane normal must be a finite non-zero vector")
        center = _as_vec3(self.center)
        if float(self.thickness) <= 0:
            raise ValueError("membrane thickness must be positive")
        object.__setattr__(self, "normal", tuple((normal / norm).tolist()))
        object.__setattr__(self, "center", tuple(center.tolist()))
        object.__setattr__(self, "thickness", float(self.thickness))
        object.__setattr__(self, "confidence", str(self.confidence).strip().lower() or "high")
        object.__setattr__(self, "method", str(self.method))
        object.__setattr__(self, "coord_frame", _normalise_frame(self.coord_frame))

    # ------------------------------------------------------------------ geometry
    @property
    def normal_array(self) -> np.ndarray:
        return np.asarray(self.normal, dtype=float)

    @property
    def center_array(self) -> np.ndarray:
        return np.asarray(self.center, dtype=float)

    @property
    def half_thickness(self) -> float:
        return 0.5 * float(self.thickness)

    def project(self, xyz: Any) -> np.ndarray:
        """Signed distance(s) of Cartesian point(s) to the mid-plane."""
        coords = np.asarray(xyz, dtype=float)
        return (coords - self.center_array) @ self.normal_array

    def compartment(self, z: Any, interface_width: float = 8.0) -> Any:
        """Classify signed distance(s) into ``core`` / ``interface`` / ``aqueous``."""
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        half = self.half_thickness
        width = float(interface_width)
        if width < 0:
            raise ValueError("interface_width must be non-negative")
        labels = np.full(z_arr.shape, "aqueous", dtype=object)
        in_core = np.abs(z_arr) < half
        in_iface = (np.abs(z_arr) >= half) & (np.abs(z_arr) < half + width)
        labels[in_core] = "core"
        labels[in_iface] = "interface"
        return labels[0] if np.isscalar(z) or np.ndim(z) == 0 else labels

    # --------------------------------------------------------------- serialise
    def as_dict(self) -> Dict[str, Any]:
        return {
            "normal": list(self.normal),
            "center": list(self.center),
            "thickness": float(self.thickness),
            "source": str(self.source),
            "confidence": str(self.confidence),
            "method": str(self.method),
            "coord_frame": str(self.coord_frame),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "MembranePlane":
        return cls(
            normal=data["normal"],
            center=data["center"],
            thickness=float(data.get("thickness", data.get("thickness_A", 30.0))),
            source=str(data.get("source", "unknown")),
            confidence=str(data.get("confidence", data.get("membrane_confidence", "high"))),
            method=str(data.get("method", data.get("membrane_method", ""))),
            coord_frame=str(data.get("coord_frame", "opm")),
        )


class MembraneMap(ABC):
    """Minimal provider interface consumed by the M sub-score."""

    name: str = "abstract"

    @abstractmethod
    def get(self, system_id: str) -> Optional[MembranePlane]:
        """Return the plane for ``system_id`` or ``None`` if unavailable."""

    def available(self, system_id: str) -> bool:
        return self.get(system_id) is not None

    def system_ids(self) -> Iterable[str]:
        return ()

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"<{type(self).__name__} name={self.name!r}>"


class NullMembraneMap(MembraneMap):
    """No membrane information: the M sub-score is neutral (1.0)."""

    name = "null"

    def get(self, system_id: str) -> Optional[MembranePlane]:
        return None


class TableMembraneMap(MembraneMap):
    """Load membrane planes from a previously exported table.

    Accepted formats
    ----------------
    CSV with columns ``system_id,nx,ny,nz,cx,cy,cz,thickness[,source]``
    JSON with ``{"planes": {system_id: {"normal": [...], "center": [...],
    "thickness": ..., "source": ...}}}``.
    """

    name = "table"

    def __init__(self, planes: Mapping[str, MembranePlane]):
        self._planes: Dict[str, MembranePlane] = {
            str(k): v for k, v in dict(planes).items()
        }

    # ------------------------------------------------------------------ provider
    def get(self, system_id: str) -> Optional[MembranePlane]:
        return self._planes.get(str(system_id))

    def system_ids(self) -> Iterable[str]:
        return tuple(self._planes)

    def __len__(self) -> int:
        return len(self._planes)

    # ------------------------------------------------------------------- loaders
    @classmethod
    def from_csv(cls, path: str | Path) -> "TableMembraneMap":
        planes: Dict[str, MembranePlane] = {}
        with Path(path).open("r", encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh):
                sid = row.get("system_id") or row.get("pdb_id") or row.get("id")
                if sid is None:
                    raise ValueError("CSV needs a system_id/pdb_id/id column")
                planes[str(sid)] = MembranePlane(
                    normal=[row["nx"], row["ny"], row["nz"]],
                    center=[row["cx"], row["cy"], row["cz"]],
                    thickness=float(row.get("thickness", 30.0) or 30.0),
                    source=str(row.get("source", "opm-table")),
                )
        return cls(planes)

    @classmethod
    def from_json(cls, path: str | Path) -> "TableMembraneMap":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        records = payload.get("planes", payload)
        return cls({str(k): MembranePlane.from_dict(v) for k, v in records.items()})

    def to_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"planes": {k: v.as_dict() for k, v in self._planes.items()}}
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return path


class OPMMembraneMap(TableMembraneMap):
    """OPM-backed plane table.

    OPM (https://opm.phar.umich.edu) has no public batch API, so the planes are
    exported once by Track C (``data/trackB/membrane_planes.json``) and loaded
    here.  ``fetch`` documents the intended extension point.
    """

    name = "opm"

    @classmethod
    def fetch(cls, system_ids: Iterable[str], cache_dir: str | Path | None = None):  # pragma: no cover
        raise NotImplementedError(
            "OPM has no public batch API. Export the PDB entries once from "
            "https://opm.phar.umich.edu and load the result with "
            "OPMMembraneMap.from_json(...) / from_csv(...)."
        )


class GeometricMembraneMap(MembraneMap):
    """Geometric fallback: build the plane from TM-helix principal axes.

    Parameters
    ----------
    systems : mapping
        ``{system_id: {"axes": (n,3) array, "midpoints": (n,3) array,
        "thickness": optional float}}`` where ``axes`` are per-helix unit
        direction vectors and ``midpoints`` are points on the helix axes.
    """

    name = "geometric"

    def __init__(self, systems: Mapping[str, Mapping[str, Any]]):
        self._systems = {str(k): dict(v) for k, v in dict(systems).items()}

    def get(self, system_id: str) -> Optional[MembranePlane]:
        data = self._systems.get(str(system_id))
        if data is None:
            return None
        axes = np.asarray(data["axes"], dtype=float)
        mids = np.asarray(data["midpoints"], dtype=float)
        if axes.ndim != 2 or axes.shape[1] != 3 or len(axes) == 0:
            raise ValueError("axes must be an (n, 3) array")
        if mids.shape != axes.shape:
            raise ValueError("midpoints must have the same shape as axes")
        axes = axes / np.linalg.norm(axes, axis=1, keepdims=True)
        # fix the arbitrary helix-axis sign by aligning every axis to the first
        ref = axes[0]
        axes = np.where((axes @ ref)[:, None] < 0, -axes, axes)
        normal = axes.mean(axis=0)
        norm = float(np.linalg.norm(normal))
        if norm < 1e-6:
            raise ValueError(f"TM-helix axes are antiparallel for {system_id!r}")
        return MembranePlane(
            normal=normal / norm,
            center=mids.mean(axis=0),
            thickness=float(data.get("thickness", 30.0)),
            source="geometric",
        )

    def system_ids(self) -> Iterable[str]:
        return tuple(self._systems)


def build_membrane_map(kind: str = "null", **kwargs: Any) -> MembraneMap:
    """Factory: ``kind`` in {null, table, opm, geometric}."""
    kind = str(kind).lower()
    if kind == "null":
        return NullMembraneMap()
    if kind in {"table", "opm"}:
        if isinstance(kwargs.get("planes"), Mapping):
            return TableMembraneMap(kwargs["planes"])
        path = kwargs.get("path")
        if path is None:
            raise ValueError("kind='table' needs path=... or planes=...")
        suffix = Path(path).suffix.lower()
        if suffix == ".csv":
            return TableMembraneMap.from_csv(path)
        return TableMembraneMap.from_json(path)
    if kind == "geometric":
        return GeometricMembraneMap(kwargs["systems"])
    raise ValueError(f"unknown membrane map kind: {kind!r}")


def load_membrane_map(path: str | Path) -> MembraneMap:
    """Auto-detect the on-disk plane table format."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        return TableMembraneMap.from_csv(path)
    return TableMembraneMap.from_json(path)


def load_c_line_planes(path: str | Path) -> TableMembraneMap:
    """Load the C-line ``membrane_planes.json`` (W2 list-of-dicts format).

    Each entry is expressed in its own OPM-oriented frame
    (``coord_frame == 'OPM_oriented'``: membrane mid-plane z = 0, normal +z), so
    planes are keyed by **pdb_id** — callers must look up by pdb_id, not by
    target.  Entries with ``source == 'unavailable'`` are skipped so that the M
    sub-score degrades to neutral (``m_neutral=True``) for those poses.
    """
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    entries = payload.get("entries", payload)
    planes: Dict[str, MembranePlane] = {}
    for entry in entries:
        source = str(entry.get("source", "c-line"))
        confidence = str(entry.get("membrane_confidence", entry.get("confidence", "high"))).strip().lower()
        if source.lower() == "unavailable" or confidence == "none":
            continue
        pid = str(entry["pdb_id"])
        planes[pid] = MembranePlane(
            normal=entry.get("normal", [0.0, 0.0, 1.0]),
            center=entry.get("center", [0.0, 0.0, 0.0]),
            thickness=float(entry.get("thickness_A", 30.0)),
            source=source,
            confidence=confidence,
            method=str(entry.get("membrane_method", "")),
            coord_frame=str(entry.get("coord_frame", "opm")),
        )
    return TableMembraneMap(planes)


def load_planes_auto(path: str | Path) -> TableMembraneMap:
    """Load either the C-line list format or the compact ``{"planes": ...}`` format."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(payload, Mapping) and "entries" in payload:
        return load_c_line_planes(path)
    return TableMembraneMap.from_json(path)


# ============================================================================
# Track B frame transforms (C-line `frame_transforms.csv`, 2026-11-16)
# ----------------------------------------------------------------------------
# Convention: x_opm = R @ x_rcsb + t, with R the 3x3 top-left block of the
# 4x4 row-major matrix (and t = column 3 / t_x,t_y,t_z).  A plane expressed in
# the OPM frame maps to the RCSB frame as:
#     normal_rcsb = R^T @ normal_opm
#     center_rcsb = R^T @ (center_opm - t)
# ============================================================================


@dataclass(frozen=True)
class FrameTransform:
    """Rigid transform from RCSB original coordinates to the OPM oriented frame."""

    pdb_id: str
    R: Any
    t: Any
    status: str = "ok"
    from_frame: str = "rcsb"
    to_frame: str = "opm"
    n_ca_matched: int = 0
    rmsd_A: float = float("nan")
    membrane_source: str = ""
    opm_plane_available: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        R = np.asarray(self.R, dtype=float).reshape(3, 3)
        t = np.asarray(self.t, dtype=float).reshape(3)
        object.__setattr__(self, "R", R)
        object.__setattr__(self, "t", t)
        object.__setattr__(self, "pdb_id", str(self.pdb_id))
        object.__setattr__(self, "status", str(self.status).strip().lower())


def load_frame_transforms(path: str | Path) -> Dict[str, FrameTransform]:
    """Load the C-line ``frame_transforms.csv`` keyed by ``pdb_id``."""
    frame = pd.read_csv(Path(path))
    out: Dict[str, FrameTransform] = {}
    for _, row in frame.iterrows():
        pid = str(row["pdb_id"])
        if str(row.get("status", "")).strip().lower() == "ok":
            R = np.array(
                [[float(row[f"r_transform_{i}{j}"]) for j in range(3)] for i in range(3)],
                dtype=float,
            )
            t = np.array(
                [float(row[f"r_transform_{i}3"]) for i in range(3)], dtype=float
            )
            if not np.all(np.isfinite(t)):
                t = np.array([float(row["t_x"]), float(row["t_y"]), float(row["t_z"])], dtype=float)
        else:  # failed transform: keep identity placeholders but flag status
            R = np.eye(3)
            t = np.zeros(3)
        out[pid] = FrameTransform(
            pdb_id=pid,
            R=R,
            t=t,
            status=str(row.get("status", "fail")),
            from_frame=str(row.get("from_frame", "rcsb")),
            to_frame=str(row.get("to_frame", "opm")),
            n_ca_matched=int(row.get("n_ca_matched", 0) or 0),
            rmsd_A=float(row.get("rmsd_A", float("nan"))),
            membrane_source=str(row.get("membrane_source", "")),
            opm_plane_available=bool(row.get("opm_plane_available", False)),
            note=str(row.get("note", "")),
        )
    return out


def transform_plane_opm_to_rcsb(plane: MembranePlane, transform: FrameTransform) -> MembranePlane:
    """Express an OPM-frame plane in the RCSB original frame (see convention above)."""
    if plane.coord_frame != "opm":
        raise ValueError(f"plane is already in frame {plane.coord_frame!r}; expected 'opm'")
    if transform.status != "ok":
        raise ValueError(f"frame transform for {transform.pdb_id} is not ok (status={transform.status})")
    R_T = transform.R.T
    normal = R_T @ plane.normal_array
    center = R_T @ (plane.center_array - transform.t)
    return MembranePlane(
        normal=normal,
        center=center,
        thickness=plane.thickness,
        source=plane.source,
        confidence=plane.confidence,
        method=plane.method,
        coord_frame="rcsb",
    )


class TrackBMembraneMap(TableMembraneMap):
    """Dual-frame membrane map for Track B (keys ``f"{pdb_id}|{frame}"``).

    * OPM-frame planes are kept as ``pdb_id|opm``;
    * when a rigid transform to RCSB is available, the same plane is also
      provided as ``pdb_id|rcsb`` (mapped with :func:`transform_plane_opm_to_rcsb`);
    * geometric planes (already RCSB_original) are exposed as ``pdb_id|rcsb``;
    * unavailable systems have no plane (M stays neutral, ``m_neutral=True``).
    ``metadata`` records per-PDB source/confidence/transform status for the
    coverage report.
    """

    name = "trackb-dual-frame"

    def __init__(self, planes: Mapping[str, MembranePlane], metadata: Mapping[str, Dict[str, Any]]):
        super().__init__(planes)
        self.metadata: Dict[str, Dict[str, Any]] = {str(k): dict(v) for k, v in dict(metadata).items()}

    def get(self, system_id: str) -> MembranePlane | None:
        """Look up ``pdb_id|frame`` first, then a frame-agnostic fallback."""
        key = str(system_id)
        if key in self._planes:
            return self._planes[key]
        # fallback: '|' separator variants / bare pdb_id for single-frame entries
        candidates = [k for k in self._planes if k.split("|")[0] == key]
        if len(candidates) == 1:
            return self._planes[candidates[0]]
        return None

    def get_for_frame(self, pdb_id: str, frame: str) -> MembranePlane | None:
        frame = _normalise_frame(frame)
        return self._planes.get(f"{pdb_id}|{frame}")

    def system_ids(self):
        return tuple(sorted({k.split("|")[0] for k in self._planes}))

    def coverage_report(self) -> Dict[str, Any]:
        """Coverage counts by membrane source / confidence / frame availability."""
        meta = self.metadata
        sources: Dict[str, int] = {}
        confidence: Dict[str, int] = {}
        transform_status: Dict[str, int] = {}
        frames: Dict[str, int] = {}
        for pid, info in meta.items():
            sources[str(info.get("membrane_source", "unknown"))] = sources.get(str(info.get("membrane_source", "unknown")), 0) + 1
            confidence[str(info.get("membrane_confidence", "none"))] = confidence.get(str(info.get("membrane_confidence", "none")), 0) + 1
            transform_status[str(info.get("transform_status", "unknown"))] = transform_status.get(str(info.get("transform_status", "unknown")), 0) + 1
            for frame in info.get("available_frames", []):
                frames[str(frame)] = frames.get(str(frame), 0) + 1
        return {
            "n_pdb_entries": len(meta),
            "by_membrane_source": sources,
            "by_confidence": confidence,
            "by_transform_status": transform_status,
            "by_available_frame": frames,
            "n_entries_with_plane": int(sum(1 for k in self._planes)),
        }


def build_trackb_membrane_map(
    planes_path: str | Path,
    transforms_path: str | Path | None = None,
) -> TrackBMembraneMap:
    """Build the dual-frame Track B membrane map (planes + frame transforms)."""
    planes = load_c_line_planes(planes_path)
    transforms = load_frame_transforms(transforms_path) if transforms_path is not None else {}
    out: Dict[str, MembranePlane] = {}
    metadata: Dict[str, Dict[str, Any]] = {}
    for pid, transform in transforms.items():
        info: Dict[str, Any] = {
            "pdb_id": pid,
            "transform_status": transform.status,
            "transform_rmsd_A": transform.rmsd_A,
            "membrane_source": transform.membrane_source or "unknown",
            "opm_plane_available": bool(transform.opm_plane_available),
            "available_frames": [],
        }
        plane = planes.get(pid)
        if plane is not None:
            info["membrane_confidence"] = str(plane.confidence)
            info["membrane_method"] = str(plane.method)
            out[f"{pid}|{plane.coord_frame}"] = plane
            info["available_frames"].append(plane.coord_frame)
            if plane.coord_frame == "opm" and transform.status == "ok":
                try:
                    out[f"{pid}|rcsb"] = transform_plane_opm_to_rcsb(plane, transform)
                    info["available_frames"].append("rcsb")
                except Exception as exc:  # pragma: no cover - defensive
                    info["transform_error"] = f"{type(exc).__name__}: {exc}"
        metadata[pid] = info
    # planes without a transform row (should not happen for Track B, but be safe)
    for pid in planes.system_ids():
        if pid in metadata or planes.get(pid) is None:
            continue
        plane = planes.get(pid)
        out[f"{pid}|{plane.coord_frame}"] = plane
        metadata[pid] = {
            "pdb_id": pid,
            "transform_status": "missing",
            "membrane_source": str(plane.source),
            "membrane_confidence": str(plane.confidence),
            "membrane_method": str(plane.method),
            "opm_plane_available": plane.coord_frame == "opm",
            "available_frames": [plane.coord_frame],
        }
    return TrackBMembraneMap(out, metadata)
