"""Part 3a: Theorem thm:main (exact count of nonreal zeros of Q_n = sum_j a_j h_{n-j}) for three admissible symbols.

For each symbol and degree n the numbers of real and nonreal zeros of Q_n are certified (certify_lib: Arb sign changes
for the real zeros, Krawczyk boxes for the nonreal ones, total = n), and compared with the number of nonreal zeros of
the symbol a in the unit disk.  The locations in (a), (b) of the theorem are measured on the approximate zeros
(z = x/sqrt n): the distance from each zeta_i = w_i + 1/w_i to the nearest zero, and how far the remaining zeros stick
out of [-2,2].  Usage: python part3a_exact_count.py [--quick]
"""
import sys

import numpy as np
from flint import arb, fmpq

from common import Report, INTERVAL, NUMERIC
from certify_lib import HermiteCombination, classify

QUICK = "--quick" in sys.argv
R = Report("part3a_exact_count")
NS = [10, 20, 40, 80, 160, 320, 400] if not QUICK else [10, 40, 160]
NS_EXTRA = {"S2": [600, 800, 1200] if not QUICK else [], "S3": [600, 800] if not QUICK else []}


def rational_series(num, pole, nmax):
    """Taylor coefficients of num(w)/(1 - w/pole) as fmpq (num: list of fmpq, pole: fmpq)."""
    q = 1 / pole
    out = []
    for j in range(nmax + 1):
        out.append(sum((num[i] * q ** (j - i) for i in range(min(j, len(num) - 1) + 1)), fmpq(0)))
    return out


def polymul(a, b):
    out = [fmpq(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


# S1: (1-2w+2w^2)(1+5w/3)/(1-10w/13): zeros (1+-i)/2 and -3/5 in D, radius 1.3
S1_num = polymul([fmpq(1), fmpq(-2), fmpq(2)], [fmpq(1), fmpq(5, 3)])
S1_zeros = [complex(0.5, 0.5), complex(0.5, -0.5), -0.6]
# S2: (1-24w/17+20w^2/17)(1+32w/13+20w^2/13)/(1-20w/21): zeros 0.6+-0.7i, -0.8+-0.1i in D, radius 1.05
S2_num = polymul([fmpq(1), fmpq(-24, 17), fmpq(20, 17)], [fmpq(1), fmpq(32, 13), fmpq(20, 13)])
S2_zeros = [complex(0.6, 0.7), complex(0.6, -0.7), complex(-0.8, 0.1), complex(-0.8, -0.1)]

symbols = {}
NMAX = max(NS + [m for v in NS_EXTRA.values() for m in v])
symbols["S1 = (1-2w+2w^2)(1+5w/3)/(1-10w/13)"] = (rational_series(S1_num, fmpq(13, 10), NMAX), S1_zeros, 1.3)
symbols["S2 = (1-24w/17+20w^2/17)(1+32w/13+20w^2/13)/(1-20w/21)"] = (rational_series(S2_num, fmpq(21, 20), NMAX), S2_zeros, 1.05)

# S3: a(w) = q(5w/6), q a sample of the Gaussian power series (degree 3000, seed 4): the symbol of eq:heat-hermite
# at t = 1.44.  Coefficients eta_j are the double-precision samples, taken as exact rationals.
eta = np.random.default_rng(4).normal(size=3001)
a3 = [fmpq(*float(e).as_integer_ratio()) * fmpq(5, 6) ** j for j, e in enumerate(eta[:NMAX + 1])]
a3f = eta[:400] * (5 / 6) ** np.arange(400)
w3 = np.roots(a3f[::-1])
w3 = w3[np.abs(w3) < 1]
# argument principle on |w|=1 and real sign changes on (-1,1) for the full degree-3000 symbol
phi = np.linspace(0, 2 * np.pi, 400001)
afull = np.polynomial.polynomial.polyval(np.exp(1j * phi), eta * (5 / 6) ** np.arange(3001))
ang = np.unwrap(np.angle(afull))
wind = int(round((ang[-1] - ang[0]) / (2 * np.pi)))
xr = np.linspace(-1 + 1e-9, 1 - 1e-9, 2000001)
ar = np.polynomial.polynomial.polyval(xr, eta * (5 / 6) ** np.arange(3001))
nreal3 = int(np.sum(np.sign(ar[1:]) != np.sign(ar[:-1])))
R.check("thm:main", "S3 (Gaussian sample at t=1.44): zeros of a in D from roots agree with the argument principle on |w|=1 "
        "and real sign changes on (-1,1)", wind == len(w3) and nreal3 == int(np.sum(np.abs(w3.imag) < 1e-12)), NUMERIC,
        f"winding number {wind}, roots in D {len(w3)}, real {nreal3}; min ||w|-1| = {np.min(np.abs(np.abs(np.roots(a3f[::-1])) - 1)):.3f}")
w3 = [complex(w.real, 0) if abs(w.imag) < 1e-12 else complex(w) for w in w3]
symbols["S3 = q(5w/6), q Gaussian sample (seed 4)"] = (a3, w3, 1.2)

table = {}
for name, (coeffs, zeros, rad) in symbols.items():
    k_nonreal = sum(1 for w in zeros if abs(w.imag) > 0)
    zetas = [w + 1 / w for w in zeros]
    rows = []
    allok = True
    extra_nonreal = {}
    for n in NS + NS_EXTRA.get(name[:2], []):
        H = HermiteCombination(lambda prec, n=n: [arb(c) for c in coeffs[:n + 1]], n)
        res = classify(H, want_full=True)
        nr = res.get("n_real_lb")
        nn = res.get("n_nonreal_lb")
        exact = res.get("exact", False)
        z = H.approx_zeros() / np.sqrt(n)
        dists = [float(np.min(np.abs(z - zt))) for zt in zetas]
        near = set()
        for zt in zetas:
            near.add(int(np.argmin(np.abs(z - zt))))
        rest = [z[i] for i in range(n) if i not in near]
        out = max([max(abs(v.real) - 2, 0) + abs(v.imag) for v in rest] + [0])
        ok_count = exact and nn == k_nonreal
        if exact and nn > k_nonreal:
            zz = z[z.imag > 1e-9]
            extra_nonreal[n] = [complex(round(v.real, 4), round(v.imag, 4)) for v in zz
                                if min(abs(v - zt) for zt in zetas) > 0.05]
        rows.append(dict(n=n, real=nr, nonreal=nn, certified=exact, predicted_nonreal=k_nonreal,
                         dist_to_zeta=[round(d, 6) for d in dists], bulk_overshoot=round(out, 6)))
        print(f"   {name[:2]} n={n}: real {nr}, nonreal {nn} (certified: {exact}); predicted nonreal {k_nonreal}; "
              f"max dist to zeta_i {max(dists):.2e}; bulk overshoot {out:.2e}", flush=True)
        allok &= bool(exact)
    table[name] = rows
    # observed onset: the smallest tested n from which the certified count equals the prediction for all larger tested n
    good = [r_["certified"] and r_["nonreal"] == k_nonreal for r_ in rows]
    n_obs = None
    for i in range(len(rows)):
        if all(good[i:]):
            n_obs = rows[i]["n"]
            break
    R.check("thm:main", f"{name}: the certified number of nonreal zeros of Q_n equals the number ({k_nonreal}) of nonreal "
            f"zeros of a in D for all tested n >= n_obs (observed onset reported; the theorem only asserts some n_0)",
            allok and n_obs is not None and good[-1] and good[-2], INTERVAL,
            f"n_obs = {n_obs}; " + "; ".join(f"n={r_['n']}: {r_['real']}R+{r_['nonreal']}C" for r_ in rows))
    if extra_nonreal:
        R.info("thm:main", f"{name}: extra nonreal zeros (upper half plane, z = x/sqrt n) before the onset",
               "; ".join(f"n={k}: {v}" for k, v in extra_nonreal.items()))
    last = rows[-1]
    R.info("thm:main(a),(b)", f"{name}: n x (distance from zeta_i to the nearest zero) and bulk overshoot, n={last['n']}",
           f"{[round(d * last['n'], 3) for d in last['dist_to_zeta']]}; overshoot {last['bulk_overshoot']}")

# independent cross-check for small n: Arb root isolation of the monomial-basis polynomial
from flint import arb_poly, acb_poly, ctx

ok = True
details = []
for name, (coeffs, zeros, rad) in list(symbols.items())[:2]:
    for n in [10, 20, 40]:
        ctx.prec = 1200
        X = arb_poly([0, 1])
        h0, h1 = arb_poly([1]), X
        hs = [h0, h1]
        for k in range(1, n):
            hs.append((X * hs[k] - arb(k).sqrt() * hs[k - 1]) * (1 / arb(k + 1).sqrt()))
        Q = sum((arb(coeffs[j]) * hs[n - j] for j in range(n + 1)), arb_poly([0]))
        rts = acb_poly(Q).roots()
        nonreal = sum(1 for r_ in rts if not r_.imag.contains(0))
        H = HermiteCombination(lambda prec, n=n: [arb(c) for c in coeffs[:n + 1]], n)
        res = classify(H, want_full=True)
        ok &= (nonreal == res.get("n_nonreal_lb")) and res.get("exact", False)
        details.append(f"{name[:2]} n={n}: {nonreal}")
R.check("thm:main", "cross-check: Arb isolation of all roots of the monomial-basis Q_n (n=10,20,40, S1 and S2) gives the "
        "same nonreal counts", ok, INTERVAL, "; ".join(details))
R.finish(extra=table)
