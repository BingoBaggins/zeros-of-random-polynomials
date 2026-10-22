"""Shared helpers for the verification suite of main/main.tex.

Every script records its checks through `Report`, which prints one line per check and writes
`results/<script>.json`.  A check has a label (the \\label in main.tex), a short description, a
method (one of SYMBOLIC, EXACT, INTERVAL, HIGHPREC-RANDOM, NUMERIC, MONTE-CARLO) and a status
(PASS / FAIL / INFO).  INFO lines carry numbers that are reported but not compared with anything.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

SYMBOLIC = "symbolic"
EXACT = "exact rational arithmetic"
INTERVAL = "interval arithmetic (Arb)"
RANDOM = "high-precision random points"
NUMERIC = "numerical (floating point)"
MC = "Monte Carlo"


class Report:
    def __init__(self, name):
        self.name = name
        self.items = []
        self.t0 = time.time()

    def check(self, label, desc, ok, method, detail=""):
        status = "PASS" if ok else "FAIL"
        self.items.append(dict(label=label, desc=desc, status=status, method=method, detail=str(detail)))
        print(f"[{status}] {label}: {desc}  ({method})" + (f"  -- {detail}" if detail else ""), flush=True)
        return ok

    def info(self, label, desc, detail, method=NUMERIC):
        self.items.append(dict(label=label, desc=desc, status="INFO", method=method, detail=str(detail)))
        print(f"[INFO] {label}: {desc}  -- {detail}", flush=True)

    def finish(self, extra=None):
        os.makedirs(RESULTS, exist_ok=True)
        nfail = sum(1 for i in self.items if i["status"] == "FAIL")
        npass = sum(1 for i in self.items if i["status"] == "PASS")
        out = dict(script=self.name, passed=npass, failed=nfail, seconds=round(time.time() - self.t0, 1),
                   checks=self.items)
        if extra is not None:
            out["data"] = extra
        with open(os.path.join(RESULTS, self.name + ".json"), "w") as f:
            json.dump(out, f, indent=1, default=str)
        print(f"== {self.name}: {npass} passed, {nfail} failed, {out['seconds']} s", flush=True)
        sys.exit(1 if nfail else 0)
