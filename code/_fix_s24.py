import csv
from pathlib import Path

rows = [
    ("method-axis amendments (frozen)",
     "plan sha256 9ce4342e45dee330243e2af64b331232b0ba381b7f5406ab41cfb30c1ca222a6; round-2 sha256 d93d8ee16bb9ff8eb7737e7e96707da1d377b3b89a40bd599dabf96a9607478b and eebedbf484bb751d651622fc2ef96b8564d785abd861032d586cd437ea3a4dff",
     "results/trackA/method_axis/FREEZE_CHAIN.json; FREEZE_CHAIN_ROUND2.json"),
    ("method-axis addenda",
     "12a: ECFP4 subset reference withdrawn (68–94% overlap); 14a: TEM1 boron B→C surrogate (22 complexes)",
     "table_S25 README; table_S28; S24 terminology"),
]
files = [Path(r"manuscript/si/table_S24_protocol_disclosure.csv"),
         Path(r"manuscript/v3/si/table_S24_protocol_disclosure.csv")]
for p in files:
    with p.open(encoding="utf-8-sig", newline="") as fh:
        data = list(csv.reader(fh))
    header = data[0]
    body = data[1:]
    body = [r for r in body if r and r[0] not in {x[0] for x in rows}]
    for r in body:
        if r and r[0] == "terminology":
            r[1] = "W12=MAES calibration gate; H1=co-folding subset; H2=predicted-receptor arm; H3=membrane MD; Addendum 12a=ECFP4 withdrawal; Addendum 14a=TEM1 surrogate"
    body.extend([list(r) for r in rows])
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(body)
    print("updated", p, "rows", len(body) + 1)
