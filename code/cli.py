"""Command line interface: ``python -m maes <command>``.

Commands
--------
smoke      run the end-to-end placeholder test (optionally write a figure)
features   compute V/M/S for a pose table (CSV in / CSV out)
fit        fit the MAES calibrator and write a model JSON
score      apply a fitted model JSON to a pose table
evaluate   independent-label evaluation of a fitted model
irt        Rasch / empirical-Bayes IRT on a (method, target, y) response table
selftest   assert-based self test (no pytest needed)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd

from . import __version__
from .api import build_features, check_monotonicity, label_audit, run_smoke
from .calibration import MAESModel
from .config import MAESConfig, default_config
from .evaluation import evaluate_predictions
from .irt import build_responses, fit_hierarchical_rasch, fit_rasch, pymc_available
from .membrane_map import NullMembraneMap, build_membrane_map, load_membrane_map
from .subscores import posebusters_available


def _dump(obj: Any, path: Optional[str]) -> None:
    text = json.dumps(obj, indent=2, default=str)
    if path:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(text, encoding="utf-8")
        print(f"[maes] wrote {path}")
    else:
        print(text)


def _load_config(path: Optional[str]) -> MAESConfig:
    return default_config() if not path else MAESConfig.load(path)


def _load_membrane(path: Optional[str], kind: str = "null"):
    if path:
        return load_membrane_map(path)
    if kind == "null":
        return NullMembraneMap()
    return build_membrane_map(kind)


def _read_table(path: str) -> pd.DataFrame:
    suffix = Path(path).suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_csv(path)


# ------------------------------------------------------------------- commands
def cmd_smoke(args: argparse.Namespace) -> int:
    report = run_smoke(
        figure_path=args.figure,
        n_poses=args.n_poses,
        n_boot=args.n_boot,
        n_weight_samples=args.n_weight_samples,
        seed=args.seed,
        verbose=True,
    )
    _dump(report, args.json_out)
    return 0


def cmd_features(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    poses = _read_table(args.poses)
    provider = _load_membrane(args.membrane_planes, kind=args.membrane_kind)
    features = build_features(poses, config=cfg, membrane_map=provider)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(out, index=False)
    print(
        f"[maes] {len(features)} rows -> {out}  "
        f"(m_neutral={features['m_neutral'].mean():.2%}, s_neutral={features['s_neutral'].mean():.2%})"
    )
    return 0


def cmd_fit(args: argparse.Namespace) -> int:
    cfg = _load_config(args.config)
    features = _read_table(args.features)
    model = MAESModel(cfg).fit(features)
    path = model.save(args.model_out)
    print(f"[maes] fitted model -> {path}")
    _dump(model.metadata, args.json_out)
    return 0


def cmd_score(args: argparse.Namespace) -> int:
    model = MAESModel.load(args.model)
    features = _read_table(args.features)
    scored = model.predict_frame(features)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(out, index=False)
    print(f"[maes] scored {len(scored)} rows -> {out}")
    return 0


def cmd_evaluate(args: argparse.Namespace) -> int:
    model = MAESModel.load(args.model)
    features = _read_table(args.features)
    metrics = evaluate_predictions(features, model)
    mono = check_monotonicity(model)
    metrics["monotonicity"] = mono
    metrics["label_audit"] = label_audit(features, model.config)
    _dump(metrics, args.json_out)
    return 0


def cmd_irt(args: argparse.Namespace) -> int:
    responses = _read_table(args.responses)
    if args.method_col and args.method_col != "method":
        responses = responses.rename(columns={args.method_col: "method"})
    if args.target_col and args.target_col != "target":
        responses = responses.rename(columns={args.target_col: "target"})
    if args.response_col and args.response_col != "y":
        responses = responses.rename(columns={args.response_col: "y"})
    result = fit_rasch(
        responses, l2=args.l2, n_boot=args.n_boot, seed=args.seed
    ) if args.model == "rasch" else fit_hierarchical_rasch(
        responses, l2=args.l2, n_boot=args.n_boot, seed=args.seed
    )
    abilities = result.ability_frame()
    difficulties = result.difficulty_frame()
    if args.out_prefix:
        prefix = Path(args.out_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        abilities.to_csv(str(prefix) + "_abilities.csv", index=False)
        difficulties.to_csv(str(prefix) + "_difficulties.csv", index=False)
        print(f"[maes] wrote {prefix}_abilities.csv / {prefix}_difficulties.csv")
    print(abilities.to_string(index=False))
    print(difficulties.to_string(index=False))
    print(f"[maes] backend={result.backend}  pymc_available={pymc_available()}  model={result.model}")
    return 0


def cmd_selftest(args: argparse.Namespace) -> int:
    from .placeholders import load_placeholders, placeholder_membrane_map

    cfg = default_config()
    issues = []

    # 1. label separation config guard
    try:
        from .config import LabelConfig, LabelSeparationError

        LabelConfig(calib_rmsd=1.0, eval_rmsd=1.0).validate()
        issues.append("label separation guard did not fire")
    except LabelSeparationError:
        pass

    # 2. end-to-end on placeholders
    poses = load_placeholders(n_poses=args.n_poses, seed=args.seed)
    features = build_features(poses, config=cfg, membrane_map=placeholder_membrane_map())
    model = MAESModel(cfg).fit(features)
    scores = model.predict(features)
    if not np.all((scores >= 0.0) & (scores <= 1.0)):
        issues.append("MAES left [0, 1]")
    mono = check_monotonicity(model)
    if not mono["ok"]:
        issues.append(f"monotonicity failed: {mono}")
    audit = label_audit(features, cfg)
    if not audit["independent"]:
        issues.append("labels are not independent")

    # 3. IRT recovers a synthetic ability ordering
    rows = []
    ability = {"m_good": 2.0, "m_mid": 0.0, "m_bad": -2.0}
    difficulty = {"t_easy": -1.5, "t_mid": 0.0, "t_hard": 1.5}
    rng = np.random.default_rng(args.seed)
    for method, theta in ability.items():
        for target, b in difficulty.items():
            p = 1.0 / (1.0 + np.exp(-(theta - b)))
            for _ in range(12):
                rows.append({"method": method, "target": target, "y": int(rng.random() < p)})
    responses = pd.DataFrame(rows)
    result = fit_rasch(responses)
    order = result.ability_frame()["method"].tolist()
    if order != ["m_good", "m_mid", "m_bad"]:
        issues.append(f"IRT ability order wrong: {order}")

    poseb = posebusters_available()
    print(f"[maes] selftest posebusters_available={poseb}  "
          f"MAES=[{scores.min():.3f},{scores.max():.3f}] rho={audit['fingerprints']['calibration'][:8]}")
    if issues:
        print("[maes] SELFTEST FAILED")
        for issue in issues:
            print("  -", issue)
        return 1
    print("[maes] SELFTEST PASS")
    return 0


# ----------------------------------------------------------------------- main
def cmd_ingest(args: argparse.Namespace) -> int:
    from .ingest import run_ingest_report

    report = run_ingest_report(
        predictions_path=args.predictions,
        manifest_path=args.manifest,
        tag=args.tag,
        write=not args.no_write,
        verbose=True,
        provisional_pb=bool(args.provisional_pb),
        pb_definition=str(args.pb_definition),
        pb_checks_path=args.pb_checks,
    )
    if report.get("report_path"):
        print(f"[maes] report -> {report['report_path']}")
    return 0 if report.get("status") in {"ok", "waiting_for_predictions", "awaiting_metrics"} else 1


# ----------------------------------------------------------------------- main
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="maes", description="Membrane-Aware Evaluation Score (v1.0)")
    parser.add_argument("--version", action="version", version=f"maes {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("smoke", help="end-to-end placeholder run")
    p.add_argument("--figure", default=None, help="write a reliability/RMSD figure here")
    p.add_argument("--json-out", default=None)
    p.add_argument("--n-poses", type=int, default=24)
    p.add_argument("--n-boot", type=int, default=0)
    p.add_argument("--n-weight-samples", type=int, default=200)
    p.add_argument("--seed", type=int, default=42)
    p.set_defaults(func=cmd_smoke)

    p = sub.add_parser("features", help="compute V/M/S from a pose table")
    p.add_argument("--poses", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--config", default=None)
    p.add_argument("--membrane-planes", default=None)
    p.add_argument("--membrane-kind", default="null")
    p.set_defaults(func=cmd_features)

    p = sub.add_parser("fit", help="fit the calibration model")
    p.add_argument("--features", required=True)
    p.add_argument("--model-out", required=True)
    p.add_argument("--config", default=None)
    p.add_argument("--json-out", default=None)
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("score", help="apply a fitted model")
    p.add_argument("--features", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_score)

    p = sub.add_parser("evaluate", help="independent-label evaluation")
    p.add_argument("--features", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--json-out", default=None)
    p.set_defaults(func=cmd_evaluate)

    p = sub.add_parser("irt", help="IRT on a response table")
    p.add_argument("--responses", required=True)
    p.add_argument("--model", choices=["rasch", "hierarchical"], default="hierarchical")
    p.add_argument("--method-col", default="method")
    p.add_argument("--target-col", default="target")
    p.add_argument("--response-col", default="y")
    p.add_argument("--l2", type=float, default=1.0)
    p.add_argument("--n-boot", type=int, default=0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-prefix", default=None)
    p.set_defaults(func=cmd_irt)

    p = sub.add_parser("ingest", help="ingest Track B predictions and print the 5-line report")
    p.add_argument("--predictions", default=None)
    p.add_argument("--manifest", default=None)
    p.add_argument("--tag", default=None)
    p.add_argument("--no-write", action="store_true")
    p.add_argument("--provisional-pb", action="store_true", help="mark strict PoseBusters bools as pending")
    p.add_argument("--pb-definition", choices=["strict", "operational"], default="strict",
                   help="PB component of the acceptable indicator (Amendment #6 operational = pb_pass_rate>=0.95)")
    p.add_argument("--pb-checks", default=None,
                   help="per-check PoseBusters matrix (default: results/trackB/pb_checks.csv if present)")
    p.set_defaults(func=cmd_ingest)

    p = sub.add_parser("selftest", help="assert-based self test")
    p.add_argument("--n-poses", type=int, default=24)
    p.add_argument("--seed", type=int, default=42)
    p.set_defaults(func=cmd_selftest)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
