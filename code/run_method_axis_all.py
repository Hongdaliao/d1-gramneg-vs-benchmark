#!/usr/bin/env python
"""Amendment #12 method-axis pipeline orchestrator (frozen protocol).

Stages (sequential, resume-safe):
  1. diffdock_main         : 806 complexes, multichain Vina-identical receptors (PO)
  2. diffdock_singlechain  : 300 complexes, GyrA/FabI/LepB chain-A single-chain sensitivity arm
  3. boltz2                : 166 nested samples (msa: empty, ligand_iptm primary)

Rules: no re-run/cherry-pick; failures recorded NA + coverage; seeds/config frozen;
DiffDock failures retried ONCE with same seed (pre-registered); logs to
results/trackA/method_axis/logs/; progress JSON updated continuously.

Run: <python> code\\trackA\\run_method_axis_all.py [--prepare-only] [--stages diffdock_main,...]
"""
import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime

ROOT = r"D:\harness\projects\d1_gramneg_vs"
BASE = os.path.join(ROOT, "results", "trackA", "method_axis")
LOGS = os.path.join(BASE, "logs")
TARGETS = ["LpxC_Ec", "GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec",
           "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
MULTICHAIN = ["GyrA_Ec", "FabI_Ec", "LepB_Ec"]
DIFFDOCK_DIR = os.path.join(ROOT, "tools", "feasibility", "diffdock", "DiffDock-main")
DIFFDOCK_PY = os.path.join(ROOT, "tools", "envs", "diffdock2", "Scripts", "python.exe")
BOLTZ_EXE = os.path.join(ROOT, "tools", "envs", "boltz", "Scripts", "boltz.exe")
BOLTZ_CACHE = os.path.join(ROOT, "tools", "boltz_work")
TORCH_HOME = os.path.join(ROOT, "tools", "feasibility", "torch_home")
DIFFDOCK_CONFIG = os.path.join(BASE, "diffdock_config.yaml")
PROGRESS = os.path.join(BASE, "progress.json")
PIDFILE = os.path.join(BASE, "pipeline.pid")
SUBSETS_DD = os.path.join(BASE, "subsets_diffdock.csv")
SUBSETS_BZ = os.path.join(BASE, "subsets_boltz2.csv")
RECEPTOR_PREP = os.path.join(BASE, "receptor_prep.json")


def now():
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def pid_alive(pid):
    try:
        import ctypes
        h = ctypes.windll.kernel32.OpenProcess(0x1000, False, int(pid))
        if not h:
            return False
        code = ctypes.c_ulong()
        ok = ctypes.windll.kernel32.GetExitCodeProcess(h, ctypes.byref(code))
        ctypes.windll.kernel32.CloseHandle(h)
        return bool(ok) and code.value == 259
    except Exception:
        return False


def ensure_dirs():
    for d in [BASE, LOGS,
              os.path.join(BASE, "diffdock", "main"),
              os.path.join(BASE, "diffdock", "singlechain"),
              os.path.join(BASE, "boltz2", "inputs"),
              os.path.join(BASE, "boltz2", "out")]:
        os.makedirs(d, exist_ok=True)


def read_csv(path):
    return list(csv.DictReader(open(path, encoding="utf-8-sig")))


def write_csv(path, rows, fields):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def write_json_atomic(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, ensure_ascii=False)
    os.replace(tmp, path)


def baseline_progress():
    if os.path.exists(PROGRESS):
        try:
            return json.load(open(PROGRESS, encoding="utf-8"))
        except Exception:
            pass
    return {"amendment": "#12 frozen 2026-11-20",
            "amendment_sha256": "b892b023f0305ece505a8925a81b8ce6beffaa8ce52809b07789a24e21ebb4e4",
            "started_at": now(), "updated_at": now(), "stage": None,
            "stages": {},
            "notes": ["protein sequence for Boltz-2 taken from receptor_prep.json (chain A) — "
                      "authoritative over subsets_boltz2.csv heuristic column"]}


def update_progress(stage=None, stage_fields=None, overall_status=None, extra=None):
    p = baseline_progress()
    if stage:
        st = p["stages"].setdefault(stage, {"status": "pending", "n_done": 0, "n_total": None,
                                            "started_at": None, "finished_at": None})
        if stage_fields:
            st.update(stage_fields)
        if st.get("started_at") is None:
            st["started_at"] = now()
    p["stage"] = stage or p.get("stage")
    p["updated_at"] = now()
    if overall_status:
        p["overall_status"] = overall_status
    if extra:
        p.update(extra)
    write_json_atomic(PROGRESS, p)
    return p


def log_line(msg, logpath=None):
    line = f"[{now()}] {msg}"
    print(line, flush=True)
    if logpath:
        with open(logpath, "a", encoding="utf-8") as f:
            f.write(line + "\n")


# ---------------------------------------------------------------- DiffDock
def make_diffdock_config():
    import yaml
    src = os.path.join(ROOT, "tools", "feasibility", "diffdock_run", "smoke_args.yaml")
    cfg = yaml.safe_load(open(src, encoding="utf-8"))
    cfg["samples_per_complex"] = 2
    cfg["inference_steps"] = 20
    cfg["actual_steps"] = 19
    import time as _t
    tmp = DIFFDOCK_CONFIG + f".tmp.{os.getpid()}"
    for attempt in range(5):
        try:
            with open(tmp, "w", encoding="utf-8") as _f:
                yaml.safe_dump(cfg, _f, sort_keys=False)
            os.replace(tmp, DIFFDOCK_CONFIG)
            return cfg
        except PermissionError:
            _t.sleep(0.5 + attempt)
    raise RuntimeError("could not write diffdock config (file locked)")


def diffdock_rows_for_target(arm, target, config):
    """Return list of input dicts for one target from the frozen subsets CSV."""
    rows = read_csv(SUBSETS_DD)
    sel = [r for r in rows if r["target"] == target]
    if arm == "singlechain" and target not in MULTICHAIN:
        return []
    if arm == "main":
        rec = None
    else:
        rec = os.path.join(BASE, "receptors_singlechain", f"{target}_singlechain.pdb")
    out = []
    for r in sel:
        protein = rec or r["receptor_pdb"]
        out.append({"target": target, "ligand_id": r["ligand_id"], "label": r["label"],
                    "smiles": r["smiles"], "protein_path": protein,
                    "complex_name": f"{target}__{r['ligand_id']}",
                    "subset_id": r["subset_id"]})
    return out


def parse_diffdock_confidence(path):
    m = re.search(r"rank1_confidence(-?\d+(?:\.\d+)?)\.sdf$", os.path.basename(path))
    return float(m.group(1)) if m else None


def parse_diffdock_output(out_dir, name):
    d = os.path.join(out_dir, name)
    if not os.path.isdir(d):
        return {"status": "failed", "confidence_rank1": None}
    files = os.listdir(d)
    rank1 = os.path.join(d, "rank1.sdf")
    conf_files = sorted([f for f in files if re.match(r"rank1_confidence-", f)])
    if not os.path.exists(rank1) or not conf_files:
        return {"status": "failed", "confidence_rank1": None}
    conf = parse_diffdock_confidence(os.path.join(d, conf_files[0]))
    if conf is None:
        return {"status": "failed", "confidence_rank1": None}
    return {"status": "ok", "confidence_rank1": conf,
            "rank1_sdf": os.path.relpath(rank1, ROOT),
            "n_rank_files": len([f for f in files if f.endswith(".sdf")])}


def run_diffdock_arm(arm, targets_filter=None, dd_limit=None):
    """Run DiffDock for one arm. Resume-safe and chunkable.

    Per target: rows with status ok/failed in scores.csv are done; pending rows are
    processed (at most dd_limit per chunk). DONE.json only when every row is
    attempted (ok or failed after the single pre-registered retry); otherwise
    PARTIAL.json records the remaining count.
    """
    stage = f"diffdock_{arm}"
    log = os.path.join(LOGS, f"{stage}.log")
    rows_all = read_csv(SUBSETS_DD)
    targets = [t for t in TARGETS if any(r["target"] == t for r in rows_all)]
    if arm == "singlechain":
        targets = [t for t in targets if t in MULTICHAIN]
    if targets_filter:
        targets = [t for t in targets if t in targets_filter]
    fields = ["arm", "target", "ligand_id", "label", "complex_name", "subset_id",
              "status", "confidence_rank1", "retry_used", "rank1_sdf", "n_rank_files"]
    # count already-processed rows across selected targets (resume-safe progress)
    n_total = sum(1 for r in rows_all if r["target"] in targets)
    n_done = 0
    for t0_ in targets:
        sp = os.path.join(BASE, "diffdock", arm, t0_, "scores.csv")
        if os.path.exists(sp):
            n_done += len(read_csv(sp))
    update_progress(stage, {"status": "running", "n_done": n_done, "n_total": n_total,
                            "chunk_limit": dd_limit})
    log_line(f"[{stage}] start targets={targets} n_total={n_total} n_done={n_done} dd_limit={dd_limit}", log)
    make_diffdock_config()
    all_rows = []
    for target in targets:
        t0 = time.time()
        tdir = os.path.join(BASE, "diffdock", arm, target)
        os.makedirs(tdir, exist_ok=True)
        rows = diffdock_rows_for_target(arm, target, None)
        if not rows:
            continue
        by_name = {r["complex_name"]: r for r in rows}
        scores_csv = os.path.join(tdir, "scores.csv")
        existing = {r["ligand_id"]: r for r in read_csv(scores_csv)} if os.path.exists(scores_csv) else {}
        # attempted = ok or failed (failed = after the pre-registered single retry)
        pending = [r for r in rows if r["ligand_id"] not in existing]
        chunk = pending[:dd_limit] if dd_limit else pending
        if not chunk:
            # everything attempted -> ensure DONE.json
            merged = [existing[r["ligand_id"]] for r in rows]
            write_csv(scores_csv, merged, fields)
            n_ok = sum(1 for x in merged if x.get("status") == "ok")
            if len(merged) == len(rows):
                write_json_atomic(os.path.join(tdir, "DONE.json"),
                                  {"arm": arm, "target": target, "n": len(merged), "n_ok": n_ok,
                                   "n_failed": len(merged) - n_ok, "elapsed_s": 0,
                                   "finished_at": now()})
            all_rows.extend(merged)
            update_progress(stage, {"n_done": sum(1 for r in all_rows)})
            log_line(f"[{stage}] {target}: nothing pending ({len(merged)}/{len(rows)} attempted)", log)
            continue
        out_dir = os.path.join(tdir, "out")
        os.makedirs(out_dir, exist_ok=True)
        # pre-scan: recover already-produced outputs from a crashed previous attempt
        pre = {}
        for r in chunk:
            res = parse_diffdock_output(out_dir, r["complex_name"])
            if res.get("status") == "ok":
                res.update({"arm": arm, "target": target, "ligand_id": r["ligand_id"],
                            "label": r["label"], "complex_name": r["complex_name"],
                            "subset_id": r["subset_id"], "retry_used": 0})
                pre[r["ligand_id"]] = res
        if pre:
            log_line(f"[{stage}] {target}: pre-scan recovered {len(pre)} existing outputs", log)
            chunk = [r for r in chunk if r["ligand_id"] not in pre]
        if not chunk:
            merged = {r["ligand_id"]: r for r in existing.values()}
            merged.update(pre)
            merged_rows = [merged[r["ligand_id"]] for r in rows if r["ligand_id"] in merged]
            write_csv(scores_csv, merged_rows, fields)
            n_ok = sum(1 for x in merged_rows if x.get("status") == "ok")
            if len(merged_rows) == len(rows):
                write_json_atomic(os.path.join(tdir, "DONE.json"),
                                  {"arm": arm, "target": target, "n": len(merged_rows), "n_ok": n_ok,
                                   "n_failed": len(merged_rows) - n_ok,
                                   "elapsed_s": round(time.time() - t0, 1), "finished_at": now()})
                log_line(f"[{stage}] {target}: TARGET DONE {n_ok}/{len(merged_rows)} ok (resume)", log)
            else:
                write_json_atomic(os.path.join(tdir, "PARTIAL.json"),
                                  {"arm": arm, "target": target, "attempted": len(merged_rows),
                                   "n_total": len(rows), "n_ok": n_ok, "updated_at": now()})
                log_line(f"[{stage}] {target}: PARTIAL {len(merged_rows)}/{len(rows)} attempted (resume)", log)
            update_progress(stage, {"n_done": len(merged_rows)})
            all_rows.extend(merged_rows)
            continue
        # build input CSV for this chunk
        inp = os.path.join(tdir, f"input_chunk_{len(existing):04d}.csv")
        write_csv(inp, [{"complex_name": r["complex_name"],
                         "protein_path": r["protein_path"],
                         "ligand_description": r["smiles"],
                         "protein_sequence": ""} for r in chunk],
                  ["complex_name", "protein_path", "ligand_description", "protein_sequence"])
        env = os.environ.copy()
        env.update({"TORCH_HOME": TORCH_HOME, "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD": "1",
                    "TEMP": os.path.join(ROOT, "tools", "tmp"), "TMP": os.path.join(ROOT, "tools", "tmp"),
                    "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
        cmd = [DIFFDOCK_PY, "-m", "inference",
               "--config", DIFFDOCK_CONFIG, "--protein_ligand_csv", inp, "--out_dir", out_dir,
               "--batch_size", "2", "--samples_per_complex", "2"]
        tlog = os.path.join(LOGS, f"{stage}_{target}.log")
        log_line(f"[{stage}] {target}: run chunk ({len(chunk)} complexes, {len(existing)} already attempted)", log)
        with open(tlog, "a", encoding="utf-8", errors="replace") as lf:
            lf.write(f"\n===== {now()} chunk {len(chunk)} =====\n")
            subprocess.run(cmd, cwd=DIFFDOCK_DIR, env=env, stdout=lf, stderr=subprocess.STDOUT,
                           timeout=6 * 3600)
        new_rows = []
        for r in chunk:
            res = parse_diffdock_output(out_dir, r["complex_name"])
            res.update({"arm": arm, "target": target, "ligand_id": r["ligand_id"],
                        "label": r["label"], "complex_name": r["complex_name"],
                        "subset_id": r["subset_id"], "retry_used": 0})
            new_rows.append(res)
        # ONE pre-registered retry (same seed) for failures in this chunk
        missing = [s for s in new_rows if s["status"] != "ok"]
        if missing:
            log_line(f"[{stage}] {target}: {len(missing)} missing in chunk -> ONE pre-registered retry", log)
            inp2 = os.path.join(tdir, f"input_chunk_{len(existing):04d}_retry.csv")
            write_csv(inp2, [{"complex_name": s["complex_name"],
                              "protein_path": by_name[s["complex_name"]]["protein_path"],
                              "ligand_description": by_name[s["complex_name"]]["smiles"],
                              "protein_sequence": ""} for s in missing],
                      ["complex_name", "protein_path", "ligand_description", "protein_sequence"])
            cmd2 = cmd.copy()
            cmd2[cmd2.index("--protein_ligand_csv") + 1] = inp2
            with open(os.path.join(LOGS, f"{stage}_{target}_retry.log"), "a",
                      encoding="utf-8", errors="replace") as lf:
                lf.write(f"\n===== {now()} retry {len(missing)} =====\n")
                subprocess.run(cmd2, cwd=DIFFDOCK_DIR, env=env, stdout=lf,
                               stderr=subprocess.STDOUT, timeout=6 * 3600)
            for s in new_rows:
                if s["status"] == "ok":
                    continue
                res = parse_diffdock_output(out_dir, s["complex_name"])
                res.update({"arm": arm, "target": target, "ligand_id": s["ligand_id"],
                            "label": s["label"], "complex_name": s["complex_name"],
                            "subset_id": s["subset_id"], "retry_used": 1})
                s.clear(); s.update(res)
        # merge and persist
        merged = {r["ligand_id"]: r for r in existing.values()}
        for s in new_rows:
            merged[s["ligand_id"]] = s
        merged_rows = [merged[r["ligand_id"]] for r in rows if r["ligand_id"] in merged]
        write_csv(scores_csv, merged_rows, fields)
        n_ok = sum(1 for x in merged_rows if x.get("status") == "ok")
        if len(merged_rows) == len(rows):
            write_json_atomic(os.path.join(tdir, "DONE.json"),
                              {"arm": arm, "target": target, "n": len(merged_rows), "n_ok": n_ok,
                               "n_failed": len(merged_rows) - n_ok,
                               "elapsed_s": round(time.time() - t0, 1), "finished_at": now()})
            log_line(f"[{stage}] {target}: TARGET DONE {n_ok}/{len(merged_rows)} ok", log)
        else:
            write_json_atomic(os.path.join(tdir, "PARTIAL.json"),
                              {"arm": arm, "target": target, "attempted": len(merged_rows),
                               "n_total": len(rows), "n_ok": n_ok, "updated_at": now()})
            log_line(f"[{stage}] {target}: CHUNK DONE {len(merged_rows)}/{len(rows)} attempted "
                     f"({n_ok} ok), {len(pending)-len(chunk)} still pending", log)
        all_rows.extend(merged_rows)
        # global progress: recount from scores files for selected targets
        done_now = 0
        for tt in targets:
            sp = os.path.join(BASE, "diffdock", arm, tt, "scores.csv")
            if os.path.exists(sp):
                done_now += len(read_csv(sp))
        update_progress(stage, {"n_done": done_now})
    # aggregate arm scores (all attempted rows for selected targets)
    agg = os.path.join(BASE, f"diffdock_{arm}_scores.csv")
    agg_rows = []
    for tt in targets:
        sp = os.path.join(BASE, "diffdock", arm, tt, "scores.csv")
        if os.path.exists(sp):
            agg_rows.extend(read_csv(sp))
    if agg_rows:
        write_csv(agg, agg_rows, fields)
    n_ok = sum(1 for r in agg_rows if r.get("status") == "ok")
    all_targets_done = all(os.path.exists(os.path.join(BASE, "diffdock", arm, tt, "DONE.json"))
                           for tt in targets)
    if all_targets_done:
        update_progress(stage, {"status": "completed", "n_done": len(agg_rows), "n_total": n_total,
                                "n_ok": n_ok, "n_failed": len(agg_rows) - n_ok,
                                "finished_at": now(), "scores_csv": os.path.relpath(agg, ROOT)})
        log_line(f"[{stage}] COMPLETE {n_ok}/{len(agg_rows)} ok -> {agg}", log)
    else:
        update_progress(stage, {"status": "running", "n_done": len(agg_rows), "n_total": n_total,
                                "n_ok": n_ok, "scores_csv": os.path.relpath(agg, ROOT)})
        log_line(f"[{stage}] CHUNK EXIT {n_ok}/{len(agg_rows)} ok (arm not complete)", log)


# ---------------------------------------------------------------- Boltz-2
def build_boltz_inputs():
    """Build one YAML per sample; protein sequence from receptor_prep (authoritative)."""
    prep = json.load(open(RECEPTOR_PREP, encoding="utf-8"))
    rows = read_csv(SUBSETS_BZ)
    out = []
    for r in rows:
        target = r["target"]
        seq = prep["targets"][target]["boltz2_sequence"]
        assert seq and len(seq) > 20, f"no sequence for {target}"
        name = f"{target}__{r['ligand_id']}"
        yp = os.path.join(BASE, "boltz2", "inputs", f"{name}.yaml")
        smi = r["smiles"].replace("'", "''")
        yaml_txt = ("version: 1\n"
                    "sequences:\n"
                    "  - protein:\n"
                    "      id: A\n"
                    f"      sequence: '{seq}'\n"
                    "      msa: empty\n"
                    "  - ligand:\n"
                    "      id: B\n"
                    f"      smiles: '{smi}'\n")
        with open(yp, "w", encoding="utf-8") as f:
            f.write(yaml_txt)
        out.append({"target": target, "ligand_id": r["ligand_id"], "label": r["label"],
                    "complex_name": name, "yaml": yp})
    meta = {"n": len(out), "built_at": now(),
            "protein_sequence_source": "receptor_prep.json boltz2_sequence (chain nearest Vina grid centre)",
            "rows": out}
    write_json_atomic(os.path.join(BASE, "boltz2", "inputs_index.json"), meta)
    return out


def parse_boltz_output(out_dir, name):
    base = os.path.join(out_dir, f"boltz_results_{name}", "predictions", name)
    if not os.path.isdir(base):
        return {"status": "failed"}
    confs = sorted([f for f in os.listdir(base) if f.startswith("confidence_") and f.endswith(".json")])
    if not confs:
        return {"status": "failed"}
    best = None
    for cf in confs:
        try:
            d = json.load(open(os.path.join(base, cf), encoding="utf-8"))
        except Exception:
            continue
        if "ligand_iptm" not in d:
            continue
        if best is None or d.get("ligand_iptm", -1) > best.get("ligand_iptm", -1):
            best = d
            best["_confidence_file"] = cf
    if best is None:
        return {"status": "failed"}
    return {"status": "ok",
            "ligand_iptm": best.get("ligand_iptm"),
            "confidence_score": best.get("confidence_score"),
            "iptm": best.get("iptm"),
            "ptm": best.get("ptm"),
            "complex_plddt": best.get("complex_plddt"),
            "complex_iplddt": best.get("complex_iplddt"),
            "pdb": os.path.relpath(os.path.join(base, best["_confidence_file"].replace("confidence_", "").replace(".json", ".pdb")), ROOT)
            if os.path.exists(os.path.join(base, best["_confidence_file"].replace("confidence_", "").replace(".json", ".pdb"))) else ""}


def run_boltz2(max_samples=None):
    stage = "boltz2"
    log = os.path.join(LOGS, "boltz2.log")
    samples = build_boltz_inputs()
    n_total = len(samples)
    update_progress(stage, {"status": "running", "n_done": 0, "n_total": n_total})
    log_line(f"[boltz2] start n={n_total}", log)
    out_root = os.path.join(BASE, "boltz2", "out")
    os.makedirs(out_root, exist_ok=True)
    scores_csv = os.path.join(BASE, "boltz2_scores.csv")
    rows = []
    if os.path.exists(scores_csv):
        prev = read_csv(scores_csv)
        rows = prev
        log_line(f"[boltz2] resume: {len(prev)} rows", log)
    # skip ok samples and samples that already exhausted the pre-registered single retry
    done_names = {r["complex_name"] for r in rows
                  if r.get("status") == "ok" or (str(r.get("attempts", "0")).isdigit() and int(r["attempts"]) >= 2)}
    n_done = len(rows)
    update_progress(stage, {"n_done": n_done})
    env = os.environ.copy()
    env.update({"BOLTZ_CACHE": BOLTZ_CACHE, "TEMP": os.path.join(ROOT, "tools", "tmp"),
                "TMP": os.path.join(ROOT, "tools", "tmp"), "PYTHONUNBUFFERED": "1",
                "PYTHONIOENCODING": "utf-8"})
    pending = [s for s in samples if s["complex_name"] not in done_names]
    if max_samples:
        pending = pending[:max_samples]
    for i, s in enumerate(pending):
        name = s["complex_name"]
        tlog = os.path.join(LOGS, f"boltz2_{name}.log")
        cmd = [BOLTZ_EXE, "predict", s["yaml"], "--out_dir", out_root,
               "--cache", BOLTZ_CACHE, "--accelerator", "gpu", "--devices", "1",
               "--output_format", "pdb", "--seed", "42", "--num_workers", "0",
               "--override", "--model", "boltz2", "--no_kernels", "--diffusion_samples", "1"]
        res = None
        for attempt in (1, 2):
            with open(tlog, "a", encoding="utf-8", errors="replace") as lf:
                lf.write(f"\n===== {now()} attempt {attempt} =====\n")
                try:
                    subprocess.run(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT, timeout=3 * 3600)
                except subprocess.TimeoutExpired:
                    lf.write("TIMEOUT\n")
            res = parse_boltz_output(out_root, name)
            if res.get("status") == "ok":
                break
        rec = {"complex_name": name, "target": s["target"], "ligand_id": s["ligand_id"],
               "label": s["label"], "status": res.get("status", "failed"),
               "ligand_iptm": res.get("ligand_iptm"), "confidence_score": res.get("confidence_score"),
               "iptm": res.get("iptm"), "ptm": res.get("ptm"),
               "complex_plddt": res.get("complex_plddt"), "complex_iplddt": res.get("complex_iplddt"),
               "pdb": res.get("pdb", ""), "attempts": attempt}
        rows = [r for r in rows if r["complex_name"] != name] + [rec]
        write_csv(scores_csv, rows, ["complex_name", "target", "ligand_id", "label", "status",
                                     "ligand_iptm", "confidence_score", "iptm", "ptm",
                                     "complex_plddt", "complex_iplddt", "pdb", "attempts"])
        n_done = len(rows)
        update_progress(stage, {"n_done": n_done})
        log_line(f"[boltz2] {i+1}/{n_total} {name}: {rec['status']} ligand_iptm={rec['ligand_iptm']}", log)
    n_ok = sum(1 for r in rows if r.get("status") == "ok")
    all_done = len(rows) >= n_total
    update_progress(stage, {"status": "completed" if all_done else "running",
                            "n_done": len(rows), "n_total": n_total,
                            "n_ok": n_ok, "n_failed": len(rows) - n_ok,
                            "finished_at": now() if all_done else None,
                            "scores_csv": os.path.relpath(scores_csv, ROOT)})
    log_line(f"[boltz2] {'COMPLETE' if all_done else 'CHUNK DONE'} {n_ok}/{len(rows)} ok -> {scores_csv}", log)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--stages", default="diffdock_main,diffdock_singlechain,boltz2")
    ap.add_argument("--dd-targets", default=None, help="comma-separated target subset for diffdock arms")
    ap.add_argument("--bz-max", type=int, default=None, help="process at most N pending Boltz-2 samples (chunk mode)")
    ap.add_argument("--dd-limit", type=int, default=None, help="process at most N pending DiffDock complexes per target (chunk mode)")
    ap.add_argument("--chunk-mode", action="store_true", help="do not mark overall run completed")
    args = ap.parse_args()
    ensure_dirs()
    make_diffdock_config()
    samples = build_boltz_inputs()
    print(f"prepared diffdock config + {len(samples)} boltz inputs")
    if args.prepare_only:
        return

    # single-instance lock (atomic create; never allow two drivers to write the same outputs)
    def acquire_lock():
        for _ in range(3):
            try:
                fd = os.open(PIDFILE, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, str(os.getpid()).encode())
                os.close(fd)
                return True
            except FileExistsError:
                try:
                    other = int(open(PIDFILE).read().strip())
                except Exception:
                    other = None
                if other and pid_alive(other):
                    print(f"another pipeline instance already running (pid={other}); exiting", flush=True)
                    return False
                try:
                    os.remove(PIDFILE)
                except OSError:
                    pass
        return False
    if not acquire_lock():
        return 0
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    try:
        for st in stages:
            dd_targets = [x.strip() for x in args.dd_targets.split(",")] if args.dd_targets else None
            if st == "diffdock_main":
                run_diffdock_arm("main", dd_targets, args.dd_limit)
            elif st == "diffdock_singlechain":
                run_diffdock_arm("singlechain", dd_targets, args.dd_limit)
            elif st == "boltz2":
                run_boltz2(args.bz_max)
            update_progress(None, None, extra={"stage": st, "overall_status": "running"})
        if args.chunk_mode:
            update_progress(None, None, overall_status="running", extra={"last_chunk_at": now()})
        else:
            update_progress(None, None, overall_status="completed", extra={"finished_at": now()})
    except Exception:
        tb = traceback.format_exc()
        with open(os.path.join(LOGS, "pipeline_error.log"), "a", encoding="utf-8") as f:
            f.write(f"\n[{now()}] ERROR\n{tb}\n")
        update_progress(None, None, overall_status="error", extra={"error": tb[-1500:]})
        raise
    finally:
        try:
            os.remove(PIDFILE)
        except OSError:
            pass


if __name__ == "__main__":
    main()
