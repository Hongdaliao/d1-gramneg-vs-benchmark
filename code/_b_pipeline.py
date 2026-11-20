#!/usr/bin/env python
"""Track A overnight pipeline: wait for LpxC_Ec, then run remaining targets in order.

Order (Lead): LpxC_Ec (running) -> GyrB_Ec -> Gyrase_Ec -> LpxC_Pa -> DHFR_Ec
              -> GyrA_Ec -> FabI_Ec -> TEM1_Ec -> LepB_Ec
Per target: manifest(20:1) -> prep(12 procs) -> prep retry + MGL fallback -> dock
            (16 workers x cpu1; exh4 for EF1 targets else exh2; box per amendment #8/#10)
            -> finalize (temporal_split/box_clipped) -> supplement rescued ligands -> finalize again
Non-EF1 (exh2): GyrA_Ec, GyrB_Ec, FabI_Ec, TEM1_Ec, LepB_Ec
Writes results/trackA/full/pipeline.log and updates target_metadata.csv / progress.csv.
Run detached: python code/trackA/_b_pipeline.py
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
LOG = os.path.join(OUTROOT, "pipeline.log")
ORDER = ["GyrB_Ec", "Gyrase_Ec", "LpxC_Pa", "DHFR_Ec", "GyrA_Ec", "FabI_Ec", "TEM1_Ec", "LepB_Ec"]
EF1 = {"LpxC_Ec", "LpxC_Pa", "DHFR_Ec", "Gyrase_Ec"}

def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def run(args, timeout=None):
    log("RUN " + " ".join(os.path.basename(a) if i == 0 else a for i, a in enumerate(args)))
    t0 = time.time()
    with open(LOG, "a", encoding="utf-8", errors="replace") as lf:
        p = subprocess.run([PY] + args, stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT, timeout=timeout)
    log(f"  -> rc={p.returncode} {time.time()-t0:.0f}s")
    return p.returncode

def wait_for(path, max_s=6 * 3600, poll=60):
    t0 = time.time()
    while time.time() - t0 < max_s:
        if os.path.exists(path):
            return True
        time.sleep(poll)
    return False

def main():
    os.makedirs(OUTROOT, exist_ok=True)
    log("pipeline start; waiting for LpxC_Ec summary.json")
    if not wait_for(os.path.join(OUTROOT, "LpxC_Ec", "summary.json")):
        log("TIMEOUT waiting for LpxC_Ec; aborting pipeline")
        return 1
    # finalize LpxC (temporal split / box_clipped / metadata) + supplement rescued ligands
    run([FIN, "--target", "LpxC_Ec"])
    run([SUP, "run", "--target", "LpxC_Ec", "--nproc", "12", "--exhaustiveness", "4"])
    run([FIN, "--target", "LpxC_Ec"])
    for target in ORDER:
        exh = "4" if target in EF1 else "2"
        log(f"=== {target} (exhaustiveness {exh}) ===")
        d = os.path.join(OUTROOT, target)
        man = os.path.join(d, "manifest.csv")
        run([V2, "manifest", "--target", target, "--ratio", "20"])
        run([V2, "prep", "--target", target, "--manifest", man, "--nproc", "12"])
        run([RETRY, "retry", "--target", target, "--nproc", "4"])
        run([MGL, "--target", target])
        run([V2, "dock", "--target", target, "--manifest", man, "--nproc", "16", "--cpu", "1", "--exhaustiveness", exh])
        run([FIN, "--target", target])
        run([SUP, "run", "--target", target, "--nproc", "12", "--exhaustiveness", exh])
        run([FIN, "--target", target])
    run([FIN, "--all"])
    log("pipeline done")
    return 0

if __name__ == "__main__":
    sys.exit(main())
