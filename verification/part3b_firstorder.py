"""Part 3b: the real-rooting time Lambda_n of P_{n,t} for fixed q, against the first-order law thm:firstorder-intro.

Lambda_n is bracketed by bisection in t (rational t): at each t the zero structure of Q_n^(t) = sum_j eta_j t^(-j/2)
h_{n-j}, which has the zeros of P_{n,t} up to scaling (eq:PQ), is certified (certify_lib): 'all-real' gives
Lambda_n <= t, a certified nonreal zero gives Lambda_n > t (lem:heat(c): the real-rooted times form [Lambda_n, oo)).
The final bracket is compared with
    n(Lambda_n - T*) ~ tau_c - W cos Psi(c_n^nat)      (thm:firstorder-intro, error O(n^-2/3), O(1/n) off the cusps)
and with the same expression at the exact phase c_n of thm:firstorder (computed from h_n(x*), h_{n-1}(x*)), and the
place of the double zero (found by Newton's method for Q = dQ/dx = 0 in (t,x), checked to lie in the bracket) with the
place formula.  Symbols: the benchmark 1+bw+2w^2 (theta*=1.1) of rem:benchmark and tab:benchmark, a rational power
series, and a sample of the Gaussian power series.  Usage: python part3b_firstorder.py [--quick]
"""
import sys

import mpmath as mp
import numpy as np
from flint import arb, fmpq, ctx

from common import Report, INTERVAL, NUMERIC
from certify_lib import HermiteCombination, classify

QUICK = "--quick" in sys.argv
R = Report("part3b_firstorder")
mp.mp.dps = 40


# ----------------------------------------------------------------------------------------------
# symbols
# ----------------------------------------------------------------------------------------------
class Symbol:
    def __init__(self, name, eta_arb, eta_mp, u_star):
        self.name = name
        self.eta_arb = eta_arb  # prec -> list of arb (length >= nmax+1)
        self.eta_mp = eta_mp  # list of mpf
        self.u = u_star  # mp complex
        self.setup()

    def q(self, w, k=0):
        # derivative of order k of the (truncated) series
        s = 0
        e = self.eta_mp
        for j in range(len(e) - 1, k - 1, -1):
            s = s * w + e[j] * mp.ff(j, k)
        return s

    def setup(self):
        u = self.u
        for _ in range(30):  # refine u_* by Newton on the full series
            u = u - self.q(u) / self.q(u, 1)
        self.u = u
        self.R = abs(u)
        self.th = mp.arg(u)
        self.T = 1 / self.R ** 2
        self.gam = u * self.q(u, 1)
        self.rho = u * self.q(u, 2) / self.q(u, 1)
        om = mp.sin(self.th)
        self.W = self.T / (2 * om ** 2)
        self.tauc = self.T / 2 * mp.cot(self.th) * mp.im(self.rho)
        self.xir = -mp.cos(self.th) * mp.re(self.rho) / 2
        self.nu = (self.th - mp.sin(self.th) * mp.cos(self.th)) / mp.pi
        self.cinf = mp.arg(self.gam) - mp.pi / 4 - self.th / 2 - mp.sin(2 * self.th) * mp.re(self.rho) / 4

    def cnat(self, n):
        return self.cinf - n * mp.pi * self.nu

    def cexact(self, n):
        xs = 2 * mp.sqrt(n) * mp.cos(self.th)
        a0, a1 = mp.mpf(1), xs
        for k in range(1, n):
            a0, a1 = a1, (xs * a1 - mp.sqrt(k) * a0) / mp.sqrt(k + 1)
        hn, hn1 = a1, a0
        car = hn + 1j * (hn * mp.cos(self.th) - hn1) / mp.sin(self.th)
        return mp.arg(car) + mp.arg(self.gam) - mp.pi / 2 - mp.sin(2 * self.th) * mp.re(self.rho) / 4


def dist_piZ(c):
    return abs(c - mp.pi * mp.nint(c / mp.pi))


def ctilde(c):
    return c - mp.pi * mp.floor(c / mp.pi + mp.mpf(1) / 2)


def Psi(c):
    d = dist_piZ(c)
    if d == 0:
        return mp.mpf(0)
    return mp.findroot(lambda P: P - mp.sin(P) - 2 * d, (mp.mpf(0), mp.pi), solver="bisect")


def F(S, c):
    return S.tauc - S.W * mp.cos(Psi(c))


def place(S, c, tau_n):
    u = S.u
    return (-mp.re(u) * mp.re(S.rho) / (2 * abs(u) ** 2) + mp.sign(ctilde(c)) * mp.sin(Psi(c)) / (2 * mp.im(u))
            + mp.re(u) * tau_n)


symbols = []
# benchmark q = 1 + b w + 2 w^2, b = -2 sqrt2 cos(1.1)
def eta_bench(prec):
    ctx.prec = prec
    return [arb(1), -2 * arb(2).sqrt() * arb("1.1").cos(), arb(2)]


bm = [mp.mpf(1), -2 * mp.sqrt(2) * mp.cos(mp.mpf("1.1")), mp.mpf(2)]
SB = Symbol("benchmark q = 1 + bw + 2w^2, theta* = 1.1", eta_bench, bm, mp.exp(1j * mp.mpf("1.1")) / mp.sqrt(2))
symbols.append(SB)

# rational q = (1+2w)(1+4w/3+20w^2/9)/(1-5w/6): u* = -0.3+0.6i, real zero -1/2, T* = 20/9
NMAX = 400
num = [fmpq(1), fmpq(4, 3) + 2, fmpq(20, 9) + fmpq(8, 3), fmpq(40, 9)]
eta_rat = [sum((num[i] * fmpq(5, 6) ** (j - i) for i in range(min(j, 3) + 1)), fmpq(0)) for j in range(1201)]
SR = Symbol("rational q = (1+2w)(1+4w/3+20w^2/9)/(1-5w/6)", lambda prec: [arb(e) for e in eta_rat],
            [mp.mpf(int(e.p)) / int(e.q) for e in eta_rat], mp.mpc(mp.mpf(-3) / 10, mp.mpf(6) / 10))
symbols.append(SR)


# Gaussian power series sample (seed chosen by the scan below), degree 4000
def gaussian_symbol(seed):
    eta = np.random.default_rng(seed).normal(size=4001)
    w = np.roots(eta[:400][::-1])
    w = w[np.abs(w) < 0.95]
    nonreal = sorted([x for x in w if x.imag > 1e-10], key=abs)
    real = [x.real for x in w if abs(x.imag) <= 1e-10]
    return eta, nonreal, real


chosen = None
for seed in range(100, 200):
    eta, nonreal, real = gaussian_symbol(seed)
    if len(nonreal) < 2:
        continue
    R0, R1 = abs(nonreal[0]), abs(nonreal[1])
    th0 = np.angle(nonreal[0])
    gap_real = min([abs(abs(x) - R0) for x in real] + [1])
    if 0.6 < R0 < 0.8 and R1 - R0 > 0.08 and gap_real > 0.08 and 0.6 < th0 < 2.5:
        chosen = (seed, eta, nonreal, real)
        break
seed, eta_g, nonreal_g, real_g = chosen
eta_g_fmpq = [fmpq(*float(e).as_integer_ratio()) for e in eta_g]
SG = Symbol(f"Gaussian power series sample (seed {seed}, degree 4000)", lambda prec: [arb(e) for e in eta_g_fmpq],
            [mp.mpf(float(e)) for e in eta_g], mp.mpc(nonreal_g[0]))
symbols.append(SG)
R.info("(G1)-(G3)", f"Gaussian sample seed {seed}: R* and the next nonreal modulus, real zeros in |w|<0.95, theta*",
       f"R*={float(SG.R):.5f}, next={abs(nonreal_g[1]):.5f}, real moduli={[round(abs(x), 4) for x in real_g]}, "
       f"theta*={float(SG.th):.4f}, T*={float(SG.T):.5f}")

# the closed forms of rem:benchmark against the generic computation
ok = (abs(SB.T - 2) < 1e-30 and abs(SB.rho - (1 - 1j * mp.cot(mp.mpf("1.1")))) < 1e-30
      and abs(SB.tauc + mp.cot(mp.mpf("1.1")) ** 2) < 1e-30 and abs(SB.tauc + SB.W - 1) < 1e-30)
R.check("rem:benchmark", "generic tau_c, W, rho* for 1+bw+2w^2 agree with the closed forms; tau_c + W = a/2 = 1", ok, NUMERIC,
        f"nu={mp.nstr(SB.nu, 6)}, tau_c={mp.nstr(SB.tauc, 6)}, W={mp.nstr(SB.W, 6)}")


# ----------------------------------------------------------------------------------------------
# certified bisection
# ----------------------------------------------------------------------------------------------
def make_H(S, n, t):
    def coeffs(prec):
        ctx.prec = prec
        e = S.eta_arb(prec)
        s = 1 / arb(t).sqrt()
        out, p = [], arb(1)
        for j in range(n + 1):
            out.append((e[j] if j < len(e) else arb(0)) * p)
            p *= s
        return out
    return HermiteCombination(coeffs, n)


def status_at(S, n, t, full=False):
    res = classify(make_H(S, n, t), want_full=full)
    if res["status"] == "undecided":
        raise RuntimeError(f"undecided at n={n}, t={t}")
    return res


def to_fmpq(x, bits=60):
    return fmpq(int(mp.nint(x * 2 ** bits)), 2 ** bits)


def bracket(S, n, tol=mp.mpf(10) ** -11):
    lo = to_fmpq(S.T + (S.tauc - S.W - 1) / n)
    hi = to_fmpq(S.T + (S.tauc + S.W + 1) / n)
    width = hi - lo
    while status_at(S, n, lo)["status"] != "nonreal":
        lo -= width
    while status_at(S, n, hi)["status"] != "all-real":
        hi += width
    steps = 0
    while float(hi - lo) > tol:
        mid = (lo + hi) / 2
        if status_at(S, n, mid)["status"] == "all-real":
            hi = mid
        else:
            lo = mid
        steps += 1
    # the counts on both sides (thm:collision(c)): exactly two nonreal zeros just before, all real just after
    before = status_at(S, n, lo, full=True)
    return lo, hi, steps, before


def newton_double_zero(S, n, t0, x0):
    """solve Q_n^(t)(x) = d_x Q_n^(t)(x) = 0 for (t,x) in mpmath."""
    e = S.eta_mp

    def QQ(t, x):
        s = 1 / mp.sqrt(t)
        a = [(e[j] if j < len(e) else 0) * s ** j for j in range(n + 1)]
        da = [-(mp.mpf(j) / 2) / t * a[j] for j in range(n + 1)]
        h = [mp.mpf(1), x]
        for k in range(1, n):
            h.append((x * h[k] - mp.sqrt(k) * h[k - 1]) / mp.sqrt(k + 1))
        Q = sum(a[j] * h[n - j] for j in range(n + 1))
        Qx = sum(a[j] * mp.sqrt(n - j) * h[n - j - 1] for j in range(n))
        Qxx = sum(a[j] * mp.sqrt((n - j) * (n - j - 1)) * h[n - j - 2] for j in range(n - 1))
        Qt = sum(da[j] * h[n - j] for j in range(n + 1))
        Qxt = sum(da[j] * mp.sqrt(n - j) * h[n - j - 1] for j in range(n))
        return Q, Qx, Qxx, Qt, Qxt

    t, x = mp.mpf(t0), mp.mpf(x0)
    for _ in range(40):
        Q, Qx, Qxx, Qt, Qxt = QQ(t, x)
        det = Qt * Qxx - Qx * Qxt
        dt = (Q * Qxx - Qx * Qx) / det
        dx = (Qt * Qx - Qxt * Q) / det
        t, x = t - dt, x - dx
        if abs(dt) < mp.mpf(10) ** -30 and abs(dx) < mp.mpf(10) ** -28:
            break
    return t, x


NS = {
    SB.name: [50, 100, 101, 110, 150, 200, 207, 250, 300, 350, 400, 408],
    SR.name: [50, 100, 150, 200, 250, 300, 350, 400],
    SG.name: [100, 150, 200, 250, 300, 350, 400],
}
EXT = [800, 1600, 3200]  # Newton's method; certified bracket [t_N - 1e-12, t_N + 1e-12] for n <= CERT_MAX
CERT_MAX = 1600
if QUICK:
    NS = {SB.name: [100, 101], SR.name: [100], SG.name: [100]}
    EXT = [800]
    CERT_MAX = 800
paper_table = {100: (0.671, 0.45390, 0.45007), 101: (0.024, -1.24112, -1.24775), 207: (1.515, 1.00171, 0.99805),
               408: (0.009, -1.37582, -1.37228)}


def predictions(S, n):
    cn, cnn = S.cexact(n), S.cnat(n)
    return cn, cnn, F(S, cn), F(S, cnn)


def double_zero(S, n, cn, pred_ex):
    y0 = 2 * mp.cos(S.th) + (S.xir + mp.sign(ctilde(cn)) * mp.sin(Psi(cn)) / (2 * mp.sin(S.th))) / n
    tN, xN = newton_double_zero(S, n, S.T + pred_ex / n, mp.sqrt(n) * y0)
    return tN, mp.sqrt(tN) * xN / mp.sqrt(n)


# in --quick mode the asymptotic checks (which assume degrees up to 3200) are reported, not judged
judge = R.check if not QUICK else (lambda label, desc, ok, method, detail="":
                                   R.info(label, desc + " [quick mode: not judged]", detail, method))

all_rows = {}
for S in symbols:
    rows = []
    for n in NS[S.name] + EXT:
        cn, cnn, pred_ex, pred_nat = predictions(S, n)
        tN, zn = double_zero(S, n, cn, pred_ex)
        if n in NS[S.name]:
            lo, hi, steps, before = bracket(S, n)
            two_before = before.get("exact") and before.get("n_nonreal_lb") == 2
            method = f"bisection ({steps} steps)"
        elif n <= CERT_MAX:
            d = fmpq(1, 10 ** 12)
            lo, hi = to_fmpq(tN, 80) - d, to_fmpq(tN, 80) + d
            b_lo = status_at(S, n, lo, full=True)
            b_hi = status_at(S, n, hi)
            if b_lo["status"] != "nonreal" or b_hi["status"] != "all-real":
                raise RuntimeError(f"certification of the Newton bracket failed at n={n}")
            two_before = b_lo.get("exact") and b_lo.get("n_nonreal_lb") == 2
            method = "Newton + certified bracket"
        else:
            lo = hi = None
            two_before = None
            method = "Newton only (not certified)"
        if lo is not None:
            lam_lo, lam_hi = mp.mpf(int(lo.p)) / int(lo.q), mp.mpf(int(hi.p)) / int(hi.q)
            tau_n = n * ((lam_lo + lam_hi) / 2 - S.T)
            in_br = bool(lam_lo <= tN <= lam_hi)
        else:
            lam_lo = lam_hi = tN
            tau_n = n * (tN - S.T)
            in_br = None
        zst = 2 * mp.cos(S.th) / S.R
        row = dict(n=n, method=method, lo=mp.nstr(lam_lo, 15), hi=mp.nstr(lam_hi, 15), two_nonreal_before=two_before,
                   dist_cn=float(dist_piZ(cn)), dist_cnat=float(dist_piZ(cnn)), n_dLambda=float(tau_n),
                   pred_exact=float(pred_ex), pred_nat=float(pred_nat), newton_in_bracket=in_br,
                   n_dz=float(n * (zn - zst)), pred_place=float(place(S, cn, tau_n)),
                   bound=float(n * (2 * mp.sqrt(mp.mpf(n) / (n - 1)) - 2)) if S is SB else None)
        rows.append(row)
        print(f"   {S.name[:9]} n={n}: Lambda in [{row['lo']}, {row['hi']}] ({method}); n(L-T*)={row['n_dLambda']:.6f}; "
              f"pred(c_n)={row['pred_exact']:.7f}, pred(c_n^nat)={row['pred_nat']:.6f}; dist={row['dist_cn']:.3f}; "
              f"n x diff={n * (row['n_dLambda'] - row['pred_nat']):.3f}; Newton in bracket {in_br}; "
              f"2 nonreal before: {two_before}", flush=True)
    all_rows[S.name] = rows
    cert = [r_ for r_ in rows if r_["newton_in_bracket"] is not None]
    ok_br = all(r_["newton_in_bracket"] and r_["two_nonreal_before"] for r_ in cert)
    R.check("thm:collision(c)", f"{S.name}: certified bracket (width <= 2e-12) for Lambda_n; exactly two nonreal zeros "
            "(certified count) at its lower end, all real at its upper end; the Newton double zero lies in the bracket",
            ok_br, INTERVAL, f"certified for n = {[r_['n'] for r_ in cert]}")
    # first-order law: scaled errors  S_n = n |diff| d_n^(1/3)  (thm:firstorder(a) says S_n = O(1))
    sc = lambda r_, key, dkey: r_["n"] * abs(r_["n_dLambda"] - r_[key]) * max(1 / r_["n"], r_[dkey]) ** (1 / 3)
    big = [r_ for r_ in rows if r_["n"] >= 800]
    s_big = max(sc(r_, "pred_nat", "dist_cnat") for r_ in big)
    s_big_ex = max(sc(r_, "pred_exact", "dist_cn") for r_ in big)
    s_all = max(sc(r_, "pred_nat", "dist_cnat") for r_ in rows)
    last = rows[-1]
    judge("thm:firstorder-intro", f"{S.name}: n(Lambda_n - T*) - (tau_c - W cos Psi(c_n^nat)) -> 0, consistent with "
            "O(n^-1 d_n^-1/3) (scaled error bounded for n >= 800, |diff| < 0.01 at n = 3200)",
            s_big < 60 and abs(last["n_dLambda"] - last["pred_nat"]) < 0.01, NUMERIC,
            f"max n|diff| d_n^(1/3): {s_big:.2f} for n>=800 (exact c_n: {s_big_ex:.2f}), {s_all:.2f} over all n; "
            f"n x diff: {[(r_['n'], round(r_['n'] * (r_['n_dLambda'] - r_['pred_nat']), 2)) for r_ in rows]}; "
            f"tau_c={float(S.tauc):.5f}, W={float(S.W):.5f}, nu={float(S.nu):.5f}")
    scp = lambda r_: r_["n"] * abs(r_["n_dz"] - r_["pred_place"]) * max(1 / r_["n"], r_["dist_cn"]) ** (2 / 3)
    judge("thm:firstorder-intro", f"{S.name}: place of the double zero, n(z_n - z*) - (formula) -> 0, consistent with "
            "O(n^-1 d_n^-2/3)", max(scp(r_) for r_ in big) < 60 and abs(last["n_dz"] - last["pred_place"]) < 0.02, NUMERIC,
            f"max n|diff| d_n^(2/3) for n>=800: {max(scp(r_) for r_ in big):.2f}; n x diff: "
            f"{[(r_['n'], round(r_['n'] * (r_['n_dz'] - r_['pred_place']), 2)) for r_ in rows]}")
    lo_end, hi_end = float(S.tauc - S.W), float(S.tauc + S.W)
    inside = all(lo_end - 0.1 < r_["n_dLambda"] < hi_end + 0.1 for r_ in big)
    judge("thm:firstorder-intro", f"{S.name}: n(Lambda_n - T*) in the limit interval [tau_c - W, tau_c + W] (up to 0.1, "
            "n >= 800)", inside, NUMERIC, f"[{lo_end:.4f}, {hi_end:.4f}], observed range over all n "
                                          f"[{min(r_['n_dLambda'] for r_ in rows):.4f}, {max(r_['n_dLambda'] for r_ in rows):.4f}]")

# the benchmark table and the elementary bound
rows = {r_["n"]: r_ for r_ in all_rows[SB.name]}
ok, det = True, []
for n, (dc, lam, pr) in paper_table.items():
    if n not in rows:
        continue
    r_ = rows[n]
    good = abs(r_["dist_cn"] - dc) < 6e-4 and abs(r_["n_dLambda"] - lam) < 6e-6 and abs(r_["pred_exact"] - pr) < 6e-6
    ok &= good
    det.append(f"n={n}: dist {r_['dist_cn']:.4f} ({dc}), n(L-2) {r_['n_dLambda']:.7f} ({lam}), "
               f"pred {r_['pred_exact']:.7f} ({pr}), n*diff {n * (r_['n_dLambda'] - r_['pred_exact']):.3f}")
R.check("tab:benchmark", "rows n=100,101,207,408 of Table tab:benchmark (dist(c_n,piZ), n(Lambda_n-T*), tau_c - W cos Psi_n "
        "with the exact c_n) reproduced to the printed digits", ok, INTERVAL, "; ".join(det))
ok = all(r_["n_dLambda"] <= r_["bound"] + 1e-9 for r_ in all_rows[SB.name] if r_["newton_in_bracket"] is not None)
r207 = rows.get(207)
R.check("rem:benchmark", "Lambda_n <= a sqrt(n/(n-1)) for every tested n; at n=207, n(Lambda_n - 2) = 1.0017 below the bound "
        "1.0036", ok and (r207 is None or (abs(r207["n_dLambda"] - 1.0017) < 1e-4)), INTERVAL,
        "; ".join(f"n={r_['n']}: {r_['n_dLambda']:.5f} <= {r_['bound']:.5f}" for r_ in all_rows[SB.name]
                  if r_["newton_in_bracket"] is not None))
vals = {n: n * abs(rows[n]["n_dLambda"] - rows[n]["pred_exact"]) for n in (110, 408) if n in rows}
if vals:
    R.check("rem:sharp", "benchmark scaled errors n|n(Lambda_n-T*) - (tau_c - W cos Psi(c_n))| near the cusps: 1.62 at n=110, "
            "1.45 at n=408", abs(vals.get(110, 1.62) - 1.62) < 0.006 and abs(vals.get(408, 1.45) - 1.45) < 0.006, INTERVAL,
            ", ".join(f"n={k}: {v:.4f} (dist(c_n,piZ) = {rows[k]['dist_cn']:.5f})" for k, v in vals.items()))
R.finish(extra=all_rows)
