#!/usr/bin/env python
"""Track A resume-all after interrupt (Lead request 2026-11-17 08:3x).

Resume semantics:
  * LpxC_Ec: rerun v2 dock -> vina_one() reuses existing poses+logs (fast), docks the rest;
    partial flush per chunk (scores.partial<chunk>.csv); then finalize + supplement + finalize.
  * Remaining targets in ORDER: skip if summary.json exists (already done), else run full steps.
  * Every target: manifest(20:1) -> prep(12) -> retry+MGL fallback -> dock(16x1, exh per EF1)
    -> finalize -> supplement -> finalize -> metadata upsert.
  * Logs to results/trackA/full/pipeline.log (timestamped).

Run detached:
  Start-Process -WindowStyle Hidden -FilePath <python> -ArgumentList code\trackA\_b_resume_all.py
"""
import json, os, subprocess, sys, time

ROOT = r"D:\harness\projects\d1_gramneg_vs"
OUTROOT = os.path.join(ROOT, "results", "trackA", "full")
PY = sys.executable
V2 = os.path.join(ROOT, "code", "trackA", "_b_trackA_v2.py")
FIN = os.path.join(ROOT, "code", "trackA", "_b_finalize.py")
SUP = os.path.join(ROOT, "code", "trackA", "_b_supplement.py")
RETRY = os.path.join(ROOT, "code", "trackA", "_b_prep_retry.py")
MGL = os.path.join(ROOT, "code", "trackA", "_b_prep_mgl.py")
LOG = os.environ.get("TRACKA_PIPELINE_LOG") or os.path.join(OUTROOT, "pipeline.log")
ORDER = ["GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
EF1 = {"LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec"}

def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    try:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except Exception:
        pass

def run(args, timeout=None):
    log("RUN " + " ".join(os.path.basename(a) if i == 0 else a for i, a in enumerate(args)))
    t0 = time.time()
    try:
        with open(LOG, "a", encoding="utf-8", errors="replace") as lf:
            p = subprocess.run([PY] + args, stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT, timeout=timeout)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        rc = -999
    log(f"  -> rc={rc} {time.time()-t0:.0f}s")
    return rc

def target_steps(target):
    d = os.path.join(OUTROOT, target)
    man = os.path.join(d, "manifest.csv")
    exh = "4" if target in EF1 else "2"
    log(f"=== {target} (exh {exh}) ===")
    run([V2, "manifest", "--target", target, "--ratio", "20"])
    run([V2, "prep", "--target", target, "--manifest", man, "--nproc", "12"])
    run([RETRY, "retry", "--target", target, "--nproc", "4"])
    run([MGL, "--target", target])
    run([V2, "dock", "--target", target, "--manifest", man, "--nproc", "16", "--cpu", "1", "--exhaustiveness", exh])
    run([FIN, "--target", target])
    run([SUP, "run", "--target", target, "--nproc", "12", "--exhaustiveness", exh])
    run([FIN, "--target", target])

LOCK = os.path.join(OUTROOT, "resume_all.lock")

def _pid_alive(pid):
    try:
        out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True, text=True, timeout=20).stdout
        return str(pid) in out
    except Exception:
        return False

def acquire_lock():
    if os.path.exists(LOCK):
        try:
            old = int(open(LOCK).read().strip() or "0")
        except Exception:
            old = 0
        if old and _pid_alive(old):
            log(f"another resume_all running (pid={old}); exit")
            return False
    open(LOCK, "w").write(str(os.getpid()))
    return True

def release_lock():
    try:
        if os.path.exists(LOCK) and int(open(LOCK).read().strip() or "0") == os.getpid():
            os.remove(LOCK)
    except Exception:
        pass

def main():
    os.makedirs(OUTROOT, exist_ok=True)
    if not acquire_lock():
        return 0
    try:
        return _main()
    finally:
        release_lock()

def _main():
    log("=== resume_all start ===")
    # 1) LpxC_Ec resume
    lx = os.path.join(OUTROOT, "LpxC_Ec", "summary.json")
    if not os.path.exists(lx):
        log("LpxC_Ec: resume dock (reuse existing poses/logs)")
        run([V2, "dock", "--target", "LpxC_Ec",
             "--manifest", os.path.join(OUTROOT, "LpxC_Ec", "manifest.csv"),
             "--nproc", "16", "--cpu", "1", "--exhaustiveness", "4"])
        run([FIN, "--target", "LpxC_Ec"])
        run([SUP, "run", "--target", "LpxC_Ec", "--nproc", "12", "--exhaustiveness", "4"])
        run([FIN, "--target", "LpxC_Ec"])
    else:
        log("LpxC_Ec already has summary.json; finalize only")
        run([FIN, "--target", "LpxC_Ec"])
    # 2) remaining targets
    for target in ORDER:
        s = os.path.join(OUTROOT, target, "summary.json")
        if os.path.exists(s):
            log(f"=== {target}: already done (summary.json) -> finalize only ===")
            run([FIN, "--target", target])
            continue
        target_steps(target)
    run([FIN, "--all"])
    log("=== resume_all done ===")

if __name__ == "__main__":
    main()
