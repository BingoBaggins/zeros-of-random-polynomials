"""Run every part of the verification suite and write results/SUMMARY.md.

    python run_all.py            # full run (about 9 minutes on a laptop)
    python run_all.py --quick    # reduced sizes (about 2 minutes)
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = ["part1_symbolic.py", "part2_constants.py", "part2_fredholm.py", "part3a_exact_count.py", "part3b_firstorder.py",
         "part4_exact_checks.py"]


def main():
    args = [a for a in sys.argv[1:] if a == "--quick"]
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    codes = {}
    for p in PARTS:
        t0 = time.time()
        print(f"=== {p}", flush=True)
        log = os.path.join(HERE, "results", p.replace(".py", ".log"))
        with open(log, "w") as f:
            proc = subprocess.run([sys.executable, os.path.join(HERE, p)] + args, cwd=HERE, stdout=f, stderr=subprocess.STDOUT)
        codes[p] = proc.returncode
        print(f"    exit {proc.returncode}, {time.time() - t0:.0f} s, log {os.path.relpath(log, HERE)}", flush=True)
    lines = ["# Verification summary", "", "| script | label | check | method | status | detail |", "|---|---|---|---|---|---|"]
    tot = {"PASS": 0, "FAIL": 0, "INFO": 0}
    for p in PARTS:
        js = os.path.join(HERE, "results", p.replace(".py", ".json"))
        if not os.path.exists(js):
            lines.append(f"| {p} | | script did not finish (exit {codes[p]}) | | FAIL | |")
            tot["FAIL"] += 1
            continue
        with open(js) as f:
            d = json.load(f)
        for c in d["checks"]:
            tot[c["status"]] += 1
            det = c["detail"].replace("|", "/")
            if len(det) > 300:
                det = det[:300] + " ..."
            lines.append(f"| {p.replace('.py', '')} | `{c['label']}` | {c['desc'].replace('|', '/')} | {c['method']} | "
                         f"{c['status']} | {det} |")
    lines[1:1] = ["", f"{tot['PASS']} passed, {tot['FAIL']} failed, {tot['INFO']} informational.", ""]
    with open(os.path.join(HERE, "results", "SUMMARY.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"{tot['PASS']} passed, {tot['FAIL']} failed, {tot['INFO']} informational; see results/SUMMARY.md")
    sys.exit(1 if tot["FAIL"] or any(codes.values()) else 0)


if __name__ == "__main__":
    main()
