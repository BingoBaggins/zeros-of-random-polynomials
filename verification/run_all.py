"""Run every part of the verification suite and write results/SUMMARY.md.

    python run_all.py            # full run (about 9 minutes on a laptop)
    python run_all.py --quick    # reduced sizes (about 2 minutes)

Requirements (tested versions): Python 3.14, sympy 1.14.0, mpmath 1.3.0, python-flint 0.9.0 (Arb ball arithmetic),
numpy 2.5.4, scipy 1.18.1.  All random inputs use fixed seeds.  The script exits nonzero if any check fails.

Contents: part1_symbolic.py (symbolic identities), part2_constants.py (the onset n_s and the constants of
prop:early, in interval arithmetic), part2_fredholm.py with fredholm_lib.py (the Fredholm determinant of
prop:fredholm), part3a_exact_count.py with certify_lib.py (certified zero counts for thm:main), part3b_firstorder.py
(Lambda_n against thm:firstorder-intro), part4_exact_checks.py (the exact computations of rem:bernoulli-general).
Each check is labelled by the \\label of the statement in the paper.  results/ holds the JSON and log output of the
last full run (140 checks passed, 0 failed).
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
