"""Verify that every R2 fill key in ``notes/maes_r2_fillkeys.md`` maps to the
formal ingest report JSON (100 % mapping check for task-25 acceptance).

Usage (from the project ``code`` directory)::

    python -m maes.fillkeys_check
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FILLKEYS = PROJECT_ROOT / "notes" / "maes_r2_fillkeys.md"
REPORT = PROJECT_ROOT / "figures" / "ingest_report_holo.json"
OUT = PROJECT_ROOT / "figures" / "fillkeys_coverage.json"

KEY_RE = re.compile(
    r"\| `((?:calibration|irt|divergence|ablation|features|exclusions)\.[^`|]+)` \|"
)


def _flatten(payload: Any, prefix: str = "") -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    if isinstance(payload, dict):
        for key, value in payload.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            flat[path] = value
            flat.update(_flatten(value, path))
    elif isinstance(payload, list):
        for idx, value in enumerate(payload):
            flat.update(_flatten(value, f"{prefix}[{idx}]"))
    return flat


def _resolve(key: str, flat: Dict[str, Any]) -> bool:
    key = key.strip()
    if key in flat:
        return True
    if key.startswith("ablation."):
        parts = key.split(".")
        variant = parts[1]
        target = f"ablation_by_variant.{variant}"
        if len(parts) == 2:
            return target in flat
        if parts[2] == "*":
            return any(path.startswith(target + ".") for path in flat)
        return f"{target}.{parts[2]}" in flat
    if key.endswith(".*"):
        return any(path.startswith(key[:-1]) for path in flat)
    return False


def check_coverage(fillkeys_path: Path = FILLKEYS, report_path: Path = REPORT) -> Dict[str, Any]:
    text = fillkeys_path.read_text(encoding="utf-8")
    keys: List[str] = []
    for line in text.splitlines():
        match = KEY_RE.search(line)
        if match:
            keys.append(match.group(1).strip())
    keys = list(dict.fromkeys(keys))
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    flat = _flatten(payload)
    mapping = {key: _resolve(key, flat) for key in keys}
    missing = [key for key, ok in mapping.items() if not ok]
    result = {
        "fillkeys": str(fillkeys_path),
        "report": str(report_path),
        "n_keys": len(keys),
        "n_mapped": int(sum(mapping.values())),
        "coverage_fraction": (float(sum(mapping.values())) / len(keys)) if keys else 1.0,
        "missing": missing,
        "mapping": mapping,
        "note": "R2 fill keys from notes/maes_r2_fillkeys.md vs figures/ingest_report_holo.json; "
                "ablation.<variant>.* resolves to ablation_by_variant.<variant>",
    }
    OUT.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    result["coverage_report"] = str(OUT)
    return result


if __name__ == "__main__":
    out = check_coverage()
    print(json.dumps({k: out[k] for k in ("n_keys", "n_mapped", "coverage_fraction", "missing", "coverage_report")}, indent=2))
