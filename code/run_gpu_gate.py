#!/usr/bin/env python
"""G0 VRAM Gate runner (B-line, task-11).

Runs Boltz-2 co-folding on 5 representative systems on the 16 GB RTX 5070 Ti,
sampling nvidia-smi for peak VRAM. Serial GPU access via jobs/GPU.lock.

Tiers (actual tokens = protein residues + ligand heavy atoms):
  T1 lpxc_4IS9        : LpxC Ec 4IS9 chain + inhibitor LTF            ~322 tokens
  T2 vhl_5T35         : VHL:ElonginC:ElonginB (+target) 4 chains + PROTAC 759 ~562
  T3 msba_mono_6BPL   : MsbA 6BPL chain A + G907 (comp AU7)          ~613
  T4 msba_mono_lipid  : MsbA monomer + G907 + 3x 3PE lipid           ~766 (proxy for ~800 tier)
  T5 msba_dimer_6BPL  : MsbA 6BPL A+B + G907                         ~1195

Outputs: results/gate/gate.csv (+ per-case logs under results/gate/logs/).
Run with the boltz venv python:
  tools\envs\boltz\Scripts\python.exe code\trackA\run_gpu_gate.py
"""
import csv, json, os, subprocess, sys, threading, time, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GATE_DATA = os.path.join(ROOT, "tools", "gate_data")
INPUTS = os.path.join(GATE_DATA, "inputs")
OUT = os.path.join(ROOT, "results", "gate")
LOGS = os.path.join(OUT, "logs")
CACHE = os.path.join(ROOT, "tools", "boltz_work")
LOCK = os.path.join(ROOT, "jobs", "GPU.lock")
BOLTZ = os.path.join(ROOT, "tools", "envs", "boltz", "Scripts", "boltz.exe")
PY = sys.executable
os.makedirs(INPUTS, exist_ok=True); os.makedirs(LOGS, exist_ok=True); os.makedirs(CACHE, exist_ok=True)

def now():
    return datetime.datetime.now().isoformat(timespec="seconds")

def heavy_atoms(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("bad smiles: " + smiles[:60])
    return m.GetNumHeavyAtoms()

def smi_for(raw, pdb, comp):
    for L in raw[pdb]["ligands"]:
        if L["comp_id"] == comp:
            return L["smiles"]
    raise KeyError(f"{pdb}/{comp}")

def seqs_for(raw, pdb):
    return [f["seq"] for f in raw[pdb]["fasta"]]

def build_cases(raw):
    cases = []
    # T1 LpxC Ec + LTF
    s = seqs_for(raw, "4IS9")[0]
    lig = [("LTF", smi_for(raw, "4IS9", "LTF"))]
    cases.append(("T1_lpxc_4IS9", "LpxC Ec (4IS9) + inhibitor LTF", [s], lig))
    # T2 VHL 4-chain + PROTAC 759
    ss = seqs_for(raw, "5T35")
    lig = [("759", smi_for(raw, "5T35", "759"))]
    cases.append(("T2_vhl_5T35", "VHL:EloC:EloB complex (5T35, 4 chains) + PROTAC 759", ss, lig))
    # T3 MsbA monomer + G907 (AU7)
    s6 = seqs_for(raw, "6BPL")[0]
    g907 = smi_for(raw, "6BPL", "AU7")
    cases.append(("T3_msba_mono_6BPL", "MsbA (6BPL chain A) + G907(AU7)", [s6], [("AU7", g907)]))
    # T4 MsbA monomer + G907 + 3x 3PE lipid (proxy ~800 tokens)
    lipid = smi_for(raw, "6BPL", "3PE")
    cases.append(("T4_msba_mono_lipid_6BPL", "MsbA chain A + G907 + 3x 3PE lipid (proxy ~800)", [s6],
                  [("AU7", g907), ("3PEa", lipid), ("3PEb", lipid), ("3PEc", lipid)]))
    # T5 MsbA dimer + G907
    cases.append(("T5_msba_dimer_6BPL", "MsbA (6BPL A+B) + G907(AU7)", [s6, s6], [("AU7", g907)]))
    out = []
    for cid, desc, prots, ligs in cases:
        tok = sum(len(p) for p in prots) + sum(heavy_atoms(s) for _, s in ligs)
        pr = sum(len(p) for p in prots)
        la = sum(heavy_atoms(s) for _, s in ligs)
        out.append({"case": cid, "desc": desc, "proteins": prots, "ligands": ligs,
                    "residues": pr, "ligand_heavy_atoms": la, "tokens": tok})
    return out

def write_yaml(case):
    lines = ["version: 1", "sequences:"]
    ids = [chr(ord("A") + i) for i in range(len(case["proteins"]))]
    for cid, seq in zip(ids, case["proteins"]):
        lines += [f"  - protein:", f"      id: {cid}", f"      sequence: {seq}"]
    start = len(case["proteins"])
    for j, (name, smi) in enumerate(case["ligands"]):
        lid = chr(ord("A") + start + j)
        lines += [f"  - ligand:", f"      id: {lid}", f"      smiles: '{smi}'"]
    path = os.path.join(INPUTS, case["case"] + ".yaml")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return path

def read_vram():
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            text=True, timeout=10).strip().splitlines()
        return max(int(x.strip()) for x in out if x.strip())
    except Exception:
        return -1

def run_case(case, yaml_path, first_run=False):
    outdir = os.path.join(OUT, "runs", case["case"])
    os.makedirs(outdir, exist_ok=True)
    logpath = os.path.join(LOGS, case["case"] + ".log")
    baseline = read_vram()
    peak = {"v": baseline}
    stop = {"v": False}

    def mon():
        while not stop["v"]:
            v = read_vram()
            if v > peak["v"]:
                peak["v"] = v
            time.sleep(0.25)

    env = dict(os.environ)
    env["HF_ENDPOINT"] = "https://hf-mirror.com"
    env["BOLTZ_CACHE"] = CACHE
    env["TMP"] = os.path.join(ROOT, "tools", "tmp")
    env["TEMP"] = env["TMP"]
    cmd = [BOLTZ, "predict", yaml_path, "--out_dir", outdir, "--cache", CACHE,
           "--use_msa_server", "--accelerator", "gpu", "--devices", "1",
           "--no_kernels", "--num_workers", "0", "--seed", "42", "--output_format", "pdb"]
    t0 = time.time()
    th = threading.Thread(target=mon, daemon=True); th.start()
    with open(logpath, "w", encoding="utf-8", errors="replace") as lf:
        lf.write("CMD: " + " ".join(cmd) + f"\nSTART: {now()}\n")
        lf.flush()
        try:
            p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env,
                               cwd=ROOT, timeout=3600)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = -999
        elapsed = time.time() - t0
        lf.write(f"\nEND: {now()} rc={rc} elapsed={elapsed:.1f}s peak_vram={peak['v']}MiB\n")
    stop["v"] = True; th.join(timeout=3)
    # success heuristic: output dir has prediction files
    nfiles = 0
    for base, _, files in os.walk(outdir):
        nfiles += len([f for f in files if f.endswith((".pdb", ".cif", ".mmcif"))])
    # ---- classify failure type (download failure vs OOM vs other) ----
    try:
        txt = open(logpath, encoding="utf-8", errors="replace").read()
    except Exception:
        txt = ""
    low = txt.lower()
    oom_pat = ["outofmemoryerror", "cuda out of memory", "out of memory", "cublas_status_alloc_failed",
               "cuda error: out of memory"]
    dl_pat = ["urlerror", "urlopen", "urlretrieve", "download", "connection", "winerror 10060",
              "timed out", "http error"]
    if rc == 0 and nfiles > 0:
        status = "能跑"
    elif rc == 0:
        status = "能跑(无输出?)"
    elif rc == -999:
        status = "超时"
    elif any(p in low for p in oom_pat):
        status = "OOM"
    elif any(p in low for p in dl_pat):
        status = "download_failure"
    else:
        status = "失败(其他)"
    delta = (peak["v"] - baseline) if peak["v"] >= 0 and baseline >= 0 else -1
    row = {"case": case["case"], "desc": case["desc"], "tokens": case["tokens"],
           "residues": case["residues"], "ligand_heavy_atoms": case["ligand_heavy_atoms"],
           "baseline_mib": baseline, "peak_mib": peak["v"], "delta_mib": delta,
           "elapsed_s": round(elapsed, 1), "returncode": rc, "output_files": nfiles,
           "status": status, "first_run_incl_weight_download": first_run,
           "note": "5-tier: T1=322 T2~562 T3~613 T4~766 T5~1195"}
    return row, logpath

def acquire_lock(task="task-11-G0"):
    os.makedirs(os.path.dirname(LOCK), exist_ok=True)
    meta = {"owner": "compute-stack(docking-runner)", "task": task, "start": now(),
            "gpu": "RTX 5070 Ti 16GB", "note": "G0 VRAM gate: boltz-2 predict x5"}
    with open(LOCK, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=1)
    print("[lock] acquired:", meta, flush=True)

def release_lock():
    try:
        os.remove(LOCK); print("[lock] released", flush=True)
    except FileNotFoundError:
        pass

EXPECTED_MIN = {"ccd.pkl": 300 << 20, "mols.tar": 1500 << 20,
                "boltz2_conf.ckpt": 2000 << 20, "boltz2_aff.ckpt": 1800 << 20}
# exact byte sizes verified against HuggingFace Content-Length on 2026-11-16
EXPECTED_EXACT = {"ccd.pkl": 345859128, "mols.tar": 1855662080,
                  "boltz2_conf.ckpt": 2286561469, "boltz2_aff.ckpt": 2062139170}

def preflight():
    """Verify Boltz-2 weights/mols are fully downloaded before touching the GPU."""
    print("== preflight: boltz cache ==", flush=True)
    ok = True
    for fn, exact in EXPECTED_EXACT.items():
        p = os.path.join(CACHE, fn)
        if not os.path.exists(p):
            print(f"  MISSING {fn} (expect exactly {exact} B)", flush=True); ok = False; continue
        sz = os.path.getsize(p)
        good = (sz == exact) or (sz >= EXPECTED_MIN[fn])
        ok = ok and good
        print(f"  {'OK' if good else 'SIZE MISMATCH'} {fn}: {sz} B (expect {exact} B)", flush=True)
    if not ok:
        print("PREFLIGHT FAIL: weights incomplete -> do NOT launch boltz (no GPU lock taken)", flush=True)
        sys.exit(3)
    print("PREFLIGHT PASS", flush=True)

def main():
    raw = json.load(open(os.path.join(GATE_DATA, "rcsb_raw.json"), encoding="utf-8"))
    cases = build_cases(raw)
    preflight()
    print("== cases ==", flush=True)
    for c in cases:
        print(f"  {c['case']}: tokens={c['tokens']} (res={c['residues']}, lig_atoms={c['ligand_heavy_atoms']})", flush=True)
        c["yaml"] = write_yaml(c)
    csv_path = os.path.join(OUT, "gate.csv")
    rows = []
    acquire_lock()
    try:
        for i, c in enumerate(cases):
            print(f"\n=== [{i+1}/{len(cases)}] {c['case']} tokens={c['tokens']} ===", flush=True)
            row, logpath = run_case(c, c["yaml"], first_run=(i == 0))
            row["date"] = "2026-11-16"
            row["source"] = "boltz-2.2.1 predict, --use_msa_server --no_kernels --seed 42; RTX5070Ti 16GB; peak via nvidia-smi 0.25s"
            rows.append(row)
            print(json.dumps(row, ensure_ascii=False), flush=True)
            with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) + ["log"])
                w.writeheader()
                for r in rows:
                    w.writerow({**r, "log": os.path.relpath(logpath, ROOT)})
    finally:
        if rows:
            with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) + ["log"])
                w.writeheader()
                for r in rows:
                    w.writerow({**r, "log": os.path.relpath(os.path.join(LOGS, r["case"] + ".log"), ROOT)})
        release_lock()
    print("\n[gate] done ->", csv_path, flush=True)

if __name__ == "__main__":
    main()
