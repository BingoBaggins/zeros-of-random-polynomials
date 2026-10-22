"""Part 2b: the law of T_* as a Fredholm determinant (prop:fredholm, rem:fredholm-num, thm:sharp, thm:lower,
rem:explicit, rem:lower, lem:qform, lem:even-bergman, lem:log-ineq).

See fredholm_lib.py for the three evaluation methods.  Usage: python part2_fredholm.py [--quick]
"""
import sys

import mpmath as mp
import numpy as np
import sympy as sp
from flint import arb, ctx, fmpq

from common import Report, SYMBOLIC, INTERVAL, NUMERIC, MC, RANDOM
from fredholm_lib import rigorous_det, varkappa_arb, nystrom_det, float_X, logdet_cond, graded_nodes

QUICK = "--quick" in sys.argv
R = Report("part2_fredholm")
rng = np.random.default_rng(20261022)
mp.mp.dps = 40
data = {}

# ----------------------------------------------------------------------------------------------
# identities used in Section 7.4
# ----------------------------------------------------------------------------------------------
err = 0
for rv in [mp.mpf("0.3"), mp.mpf("0.77"), mp.mpf("0.95")]:
    x = rv * rv
    S = {k: mp.quad(lambda t: mp.sin(k * t) / abs(1 - x * mp.exp(2j * t)), [0, mp.pi / 2, mp.pi]) for k in range(-1, 26, 2)}
    err = max(err, abs(S[1] - 2 * mp.atan(rv) / rv), abs(S[-1] + S[1]))
    for k in range(1, 22, 2):
        err = max(err, abs((k + 1) * x * S[k + 2] - k * (1 + x ** 2) * S[k] + (k - 1) * x * S[k - 2] + 2 * (1 - x)))
    c = [mp.mpf(1)]
    Nx = int(mp.ceil(mp.log(mp.mpf(10) ** -34) / mp.log(x)))
    while len(c) <= Nx:
        c.append(c[-1] * (2 * len(c) - 1) / (2 * len(c)))
    for k in [1, 5, 13]:
        ser = sum(c[a] * c[s_ - a] * x ** s_ * 2 * k / (k * k - 4 * (2 * a - s_) ** 2) for s_ in range(Nx + 1) for a in range(s_ + 1))
        err = max(err, abs(ser - S[k]))
R.check("rem:fredholm-num", "S_1 = 2 arctan(rho)/rho, S_-k = -S_k, the three-term relation for S_k (k<=21), and the series "
        "S_k = sum c_a c_b x^(a+b) 2k/(k^2-4(a-b)^2) used below, at rho in {0.3,0.77,0.95}", err < mp.mpf(10) ** -25, RANDOM,
        f"max error {mp.nstr(err, 3)}")
xs = sp.symbols("x")
ser1 = sp.series(2 * xs ** 2 / (1 - xs ** 2) ** 2, xs, 0, 40).removeO()
ser2 = sp.series((1 + xs ** 2) / (1 - xs ** 2) ** 2, xs, 0, 40).removeO()
ok = sp.expand(ser1 - sum(m * xs ** m for m in range(0, 40, 2))) == 0 and sp.expand(ser2 - sum((m + 1) * xs ** m for m in range(0, 40, 2))) == 0
ser3 = sp.series(2 * xs / (1 - xs ** 2) ** 2, xs, 0, 40).removeO()
ok &= sp.expand(ser3 - sum(2 * m * xs ** (2 * m - 1) for m in range(1, 21))) == 0  # sum |A_mn| r^(m+n)
vs = sp.symbols("varsigma", positive=True)
ok &= sp.simplify((sp.exp(-vs) / (1 - sp.exp(-2 * vs)) - 1 / (2 * sp.sinh(vs))).rewrite(sp.exp)) == 0  # geometric series
R.check("eq:Lkernel", "sum_{m even} m x^m = 2x^2/(1-x^2)^2, sum (m+1)x^m = (1+x^2)/(1-x^2)^2, sum|A_mn|r^(m+n) = 2r/(1-r^2)^2, "
        "eq:traces sum_k e^{-(2k+1)s} = 1/(2 sinh s)", ok, SYMBOLIC)
v = sp.symbols("v")
ok = sp.simplify(sp.integrate((1 - sp.sin(v)) ** 2 / sp.cos(v) ** 2, (v, 0, sp.pi / 2)) - (2 - sp.pi / 2)) == 0
ok &= abs(mp.quad(lambda u: -mp.log(1 - mp.exp(-u)), [0, 1, mp.inf]) - mp.pi ** 2 / 6) < mp.mpf(10) ** -30
R.check("lem:even-bergman / thm:lower", "int_0^{pi/2} (1-sin v)^2/cos^2 v dv = 2-pi/2 (symbolic); int_0^oo -log(1-e^-u) du = "
        "pi^2/6 (40-digit quadrature)", ok, SYMBOLIC)
err = 0
for _ in range(20):
    u_, v_ = rng.uniform(-3, 3), rng.uniform(0.01, np.pi / 2 - 0.01)
    zeta = u_ + 1j * v_
    w = np.tanh(zeta / 2)
    dw = abs(0.5 / np.cosh(zeta / 2) ** 2)
    den = np.cosh(u_) + np.cos(v_)
    tau = 2 * w.imag / abs(1 - w * w)
    err = max(err, abs(dw - 1 / den), abs(abs(1 - w * w) - 2 / den), abs(1 - abs(w) ** 2 - 2 * np.cos(v_) / den),
              abs(tau - np.sin(v_)), max(0, np.tanh(u_ / 2) ** 2 - abs(w) ** 2))
R.check("lem:even-bergman", "strip substitution w = tanh(zeta/2): |dw/dzeta|, |1-w^2|, 1-|w|^2, tau = sin v, "
        "|w|^2 >= tanh^2(u/2)", err < 1e-12, NUMERIC, f"max error {err:.2e}")
worst = 0
for p in range(2, 9):
    Cp = max(2 * np.e, 1 + sum(1 / j for j in range(1, p)))
    lam = (rng.uniform(0, 3, 20000) * np.exp(1j * rng.uniform(0, 2 * np.pi, 20000)))
    lam = np.concatenate([lam, 1 + 1e-3 * np.exp(1j * rng.uniform(0, 2 * np.pi, 2000))])
    lhs = np.log(np.abs(1 - lam))
    rhs = -np.real(sum(lam ** j / j for j in range(1, p))) + Cp * np.abs(lam) ** p
    worst = max(worst, np.max(lhs - rhs))
R.check("lem:log-ineq", "log|1-l| <= -Re sum_{j<p} l^j/j + C_p|l|^p, p=2..8, 154000 random l", worst <= 1e-12, NUMERIC,
        f"max(lhs-rhs) = {worst:.3e}")

# ----------------------------------------------------------------------------------------------
# (1) rigorous enclosures
# ----------------------------------------------------------------------------------------------
cases = [((1, 2), 60, 80), ((7, 10), 110, 140), ((4, 5), 170, 210), ((9, 10), 400, 430)]
if QUICK:
    cases = cases[:2]
rig = {}
for rf, M, N in cases:
    res = rigorous_det(rf, M, N)
    r = arb(fmpq(*rf))
    rig[rf] = res
    d = res["det"]
    lohi = f"[{d.lower().str(14, radius=False)}, {d.upper().str(14, radius=False)}]"
    data[f"rigorous r={rf[0]}/{rf[1]}"] = dict(det=lohi, M=M, N=N, tail=res["tauM"].str(3), invnorm=res["invnorm"].str(5),
                                               trace=res["trace"].str(20))
    R.info("prop:fredholm", f"r={rf[0]}/{rf[1]} (t={float(1 / r ** 2):.6f}): P(T*<=t) = det(I-X(r))",
           f"{lohi}  [M={M}, N={N}, truncation tail {res['tauM'].mid().str(3)}, ||(I-X_M)^-1||<={res['invnorm'].mid().str(4)}]",
           INTERVAL)
    R.check("prop:fredholm", f"tr X(r) = varkappa(r) at r={rf[0]}/{rf[1]}", res["trace"].overlaps(varkappa_arb(r)), INTERVAL,
            f"tr X = {res['trace'].str(15)}, varkappa = {varkappa_arb(r).str(15)}")
    prod = arb(1)
    for k in range(0, 4000):
        prod *= 1 - r ** (4 * k + 1)
    R.check("thm:lower", f"P(T*<=t) >= prod_k (1-t^(-2k-1/2)) at r={rf[0]}/{rf[1]}", d > prod, INTERVAL,
            f"{d.str(8)} > {prod.str(8)}")
if (9, 10) in rig:
    d = rig[(9, 10)]["det"]
    claim = arb("0.10868062")
    ok = (d - claim).abs_upper() < arb("5e-9")
    R.check("rem:fredholm-num", "P(T* <= 100/81) = 0.10868062 (r=0.9)", ok, INTERVAL,
            f"certified enclosure [{d.lower().str(12, radius=False)}, {d.upper().str(12, radius=False)}]")
    prod = arb(1)
    for k in range(0, 2000):
        prod *= 1 - arb("0.9") ** (4 * k + 1)
    R.check("rem:lower(b)", "the bound of thm:lower at t=100/81 is 0.0112", (prod - arb("0.0112")).abs_upper() < arb("0.00005"),
            INTERVAL, f"prod = {prod.str(6)}")

# ----------------------------------------------------------------------------------------------
# (2) Bornemann quadrature (Nystrom) on the half disk
# ----------------------------------------------------------------------------------------------
for rf, grids in [((1, 2), [(8, 16), (12, 24), (16, 32), (24, 48)]),
                  ((9, 10), [(16, 32), (24, 48), (32, 64), (40, 80), (48, 96)] if not QUICK else [(16, 32), (24, 48)])]:
    rv = rf[0] / rf[1]
    ref = float(rig[rf]["det"].mid()) if rf in rig else None
    rows = []
    for nr, nt in grids:
        dval, ld, sg = nystrom_det(rv, nr, nt)
        rows.append((nr, nt, dval.real, abs(dval.imag), (dval.real - ref) if ref is not None else None))
    data[f"nystrom r={rv}"] = rows
    best = rows[-1]
    R.info("rem:fredholm-num", f"Bornemann/Nystrom, half disk, r={rv}: (n_rho, n_theta, det, |Im det|, det - certified value)",
           "; ".join(f"({a},{b}): {c:.12f}, {d_:.1e}, {e:+.1e}" for a, b, c, d_, e in rows), NUMERIC)
    if ref is not None:
        R.check("rem:fredholm-num", f"Bornemann quadrature reproduces the certified value at r={rv}", abs(best[4]) < 1e-9,
                NUMERIC, f"difference {best[4]:.2e} with {best[0]}x{best[1]} nodes")

# ----------------------------------------------------------------------------------------------
# (3) double-precision matrix method, up to r = 0.99
# ----------------------------------------------------------------------------------------------
flt = {}
for rv in ([0.9, 0.95, 0.99] if not QUICK else [0.9]):
    M0 = int(14 / (1 - rv))
    M0 += M0 % 2
    out = []
    for M in [M0, int(1.25 * M0) // 2 * 2]:
        X = float_X(rv, M)
        sg, ld, kap = logdet_cond(X)
        out.append((M, sg, ld, kap, np.trace(X)))
    flt[rv] = (out, X)
    (M1, s1, l1, k1, tr1), (M2, s2, l2, k2, tr2) = out
    kr = float(varkappa_arb(arb(rv)).mid()) if rv != 0.99 else float(varkappa_arb(arb("0.99")).mid())
    data[f"float r={rv}"] = dict(M=[M1, M2], minus_log_P=[-l1, -l2], kappa_det=[k1, k2], trace=[tr1, tr2], varkappa=kr)
    R.info("rem:fredholm-num", f"double precision, r={rv}: -log P(T*<=t) for M={M1},{M2}; componentwise condition number "
           "of det; tr X - varkappa", f"{-l1:.9f}, {-l2:.9f}; kappa={k1:.1f}; trace error {tr2 - kr:.1e}", NUMERIC)
    if rv == 0.9 and (9, 10) in rig:
        ref = -float(rig[(9, 10)]["det"].mid().log().mid()) if False else -np.log(float(rig[(9, 10)]["det"].mid()))
        R.check("rem:fredholm-num", "double-precision matrix method agrees with the certified value at r=0.9", abs(-l2 - ref) < 1e-9,
                NUMERIC, f"-log P: {-l2:.12f} vs certified {ref:.12f}")
    if rv == 0.99:
        R.check("rem:fredholm-num", "-log P(T* <= 10000/9801) = 36.180955 (r=0.99)", abs(-l2 - 36.180955) < 2e-6 and abs(l1 - l2) < 1e-6,
                NUMERIC, f"{-l1:.8f} (M={M1}), {-l2:.8f} (M={M2}); det condition number {k2:.0f}, so an entrywise relative "
                         f"error e changes -log P by at most ~{k2:.0f} e")
        S = (X + X.T) / 2
        sgS, ldS = np.linalg.slogdet(np.eye(len(S)) - S)
        lprod = float(sum(np.log1p(-0.99 ** (4 * k + 1)) for k in range(20000)))
        R.check("rem:lower(b)", "at r=0.99: log P = -36.181, log det(I-S) = -36.269, sum log(1-r^(4k+1)) = -42.090",
                abs(l2 + 36.181) < 5e-4 and abs(ldS + 36.269) < 5e-4 and abs(lprod + 42.090) < 5e-4 and sgS > 0, NUMERIC,
                f"{l2:.4f}, {ldS:.4f}, {lprod:.4f}")
        ev = np.linalg.eigvals(X)
        R.info("q:bernoulli", "eigenvalues of X(0.99) (M=%d): max |Im|, min/max Re of those with |l|>1e-12" % M2,
               f"{np.max(np.abs(ev.imag)):.2e}, {ev.real[np.abs(ev) > 1e-12].min():.3e}, {ev.real.max():.6f}")
        tt = 1 / 0.99 ** 2
        R.info("thm:sharp", "(t-1)(-log P(T*<=t)) at r=0.99, against pi^2/12 = 0.822467", f"{(tt - 1) * (-l2):.6f}")

# ----------------------------------------------------------------------------------------------
# (4) quadratic form, even Bergman kernel, explicit bound
# ----------------------------------------------------------------------------------------------
if (9, 10) in rig:
    rv = 0.9
    Xr = rig[(9, 10)]["X"]
    nM = Xr.nrows()
    Xf = np.array([[float(Xr[i, j].mid()) for j in range(nM)] for i in range(nM)])
    xg, wg = np.polynomial.legendre.leggauss(400)
    th = np.pi * (xg + 1) / 2
    wt = np.pi / 2 * wg
    xr_ = rv * xg
    wr_ = rv * wg
    sig = (1 + rv ** 2) * np.sin(th) / np.abs(1 - rv ** 2 * np.exp(2j * th))
    worst, smax = 0, sig.max()
    for _ in range(20):
        cvec = np.zeros(nM)
        cvec[:11] = rng.normal(size=11)
        lhs = cvec @ Xf @ cvec
        fz = sum(cvec[i] * (rv * np.exp(1j * th)) ** (2 * i) for i in range(11))
        fx = sum(cvec[i] * xr_ ** (2 * i) for i in range(11))
        rhs = rv / np.pi * np.sum(wt * np.abs(fz) ** 2 * sig) - 1 / np.pi * np.sum(wr_ * fx ** 2)
        worst = max(worst, abs(lhs - rhs) / (abs(lhs) + 1e-300))
    R.check("lem:qform", "c^T X(r) c = (r/pi) int |f(re^it)|^2 sigma_r - (1/pi) int_{-r}^r f^2, r=0.9, 20 random c "
            "(m<=20); 0 <= sigma_r <= 1", worst < 1e-10 and smax <= 1 + 1e-12, NUMERIC, f"max rel. difference {worst:.1e}")
    Sn = (Xf + Xf.T) / 2
    Rn = np.diag(rv ** (2 * np.arange(0, 2 * nM, 2) + 1))
    R.check("thm:lower", "S_n <= R_n (Step 1) for the compression of X(0.9) to m<=%d" % (2 * nM - 2),
            np.linalg.eigvalsh(Rn - Sn).min() > -1e-12, NUMERIC, f"min eigenvalue of R_n - S_n = {np.linalg.eigvalsh(Rn - Sn).min():.2e}")


def even_bergman_norm(rv, n1=40, n2=80):
    xg, wg = np.polynomial.legendre.leggauss(n1)
    rho = rv * (xg + 1) / 2
    wr = rv / 2 * wg
    xt, wt_ = np.polynomial.legendre.leggauss(n2)
    thf = np.pi * (xt + 1)  # full circle for z
    wtf = np.pi * wt_
    thh = np.pi * (xt + 1) / 2  # upper half for w
    wth = np.pi / 2 * wt_
    Z = (rho[:, None] * np.exp(1j * thf[None, :])).ravel()
    WZ = (wr[:, None] * rho[:, None] * wtf[None, :]).ravel()
    Wp = (rho[:, None] * np.exp(1j * thh[None, :])).ravel()
    WW = (wr[:, None] * rho[:, None] * wth[None, :]).ravel()
    tot = 0.0
    for i0 in range(0, len(Wp), 400):
        wv = Wp[i0:i0 + 400][None, :]
        wb = np.conj(wv)
        zc = Z[:, None]
        Lk = 1 / (2j * np.pi * np.abs(1 - wv ** 2)) * ((zc - wb) / (1 - zc * wb) ** 2 - (zc + wb) / (1 + zc * wb) ** 2)
        Kk = 1 / (2 * np.pi) * (1 / (1 - zc * wb) ** 2 + 1 / (1 + zc * wb) ** 2)
        tot += np.sum(WZ[:, None] * WW[i0:i0 + 400][None, :] * np.abs(Lk - Kk) ** 2)
    return np.sqrt(2 * tot)


for rv in [0.5, 0.8, 0.9]:
    nrm = even_bergman_norm(rv)
    nrm2 = even_bergman_norm(rv, 50, 100)
    bound = 1 + np.log(1 / (1 - rv))
    R.check("lem:even-bergman", f"||L - K||_2 <= 1 + log(1/(1-r)) at r={rv}", nrm2 <= bound, NUMERIC,
            f"||L-K||_2 = {nrm2:.6f} (coarser grid {nrm:.6f}), bound {bound:.4f}")

# rem:explicit at r = 0.8: (t-1)(varkappa - ||L||_2^2/2) = 0.2796, and exp(-varkappa+||L||^2/2) >= P(T*<=t)
rv = 0.8
xg, wg = np.polynomial.legendre.leggauss(120)
rho = rv * (xg + 1) / 2
wr = rv / 2 * wg
xt, wt_ = np.polynomial.legendre.leggauss(240)
thf = np.pi * (xt + 1)
wtf = np.pi * wt_
Z = (rho[:, None] * np.exp(1j * thf[None, :])).ravel()
WZ = (wr[:, None] * rho[:, None] * wtf[None, :]).ravel()
L2 = 0.0
for m in range(0, 400, 2):
    psi = m * Z ** (m - 1) - (m + 1) * Z ** (m + 1) if m > 0 else -(m + 1) * Z ** (m + 1)
    vm2 = np.sum(WZ * np.abs(psi) ** 2 / (np.pi ** 2 * np.abs(1 - Z ** 2) ** 2))
    L2 += np.pi * rv ** (2 * m + 2) / (m + 1) * vm2
kap8 = float(varkappa_arb(arb("0.8")).mid())
tt = 1 / rv ** 2
val = (tt - 1) * (kap8 - L2 / 2)
ok = abs(val - 0.2796) < 5e-5
if (4, 5) in rig:
    ok &= np.exp(-kap8 + L2 / 2) >= float(rig[(4, 5)]["det"].mid())
R.check("rem:explicit", "(t-1)(varkappa(r) - ||L||_2^2/2) = 0.2796 at r=0.8, and P(T*<=t) <= exp(-varkappa + ||L||_2^2/2)", ok,
        NUMERIC, f"value {val:.5f}; exp bound {np.exp(-kap8 + L2 / 2):.6f} vs P = "
                 f"{float(rig[(4, 5)]['det'].mid()) if (4, 5) in rig else float('nan'):.6f}")

# ----------------------------------------------------------------------------------------------
# (5) Monte Carlo: the law of k_r (number of zeros of q in the upper half disk)
# ----------------------------------------------------------------------------------------------
mc_cases = [((1, 2), 45, 20000), ((7, 10), 80, 12000), ((9, 10), 230, 6000)]
if QUICK:
    mc_cases = [((1, 2), 45, 3000)]
for rf, D, nsamp in mc_cases:
    rv = rf[0] / rf[1]
    counts = np.zeros(nsamp, dtype=int)
    for i in range(nsamp):
        eta = rng.normal(size=D + 1)
        z = np.roots(eta[::-1])
        counts[i] = np.sum((np.abs(z) < rv) & (z.imag > 0))
    p0, ek = np.mean(counts == 0), counts.mean()
    if rf in rig:
        Xr = rig[rf]["X"]
        nM = Xr.nrows()
        Xf = np.array([[float(Xr[i, j].mid()) for j in range(nM)] for i in range(nM)])
        lam = np.linalg.eigvals(Xf)
        # E s^k = prod (1 - (1-s) l): coefficients in s
        poly = np.array([1.0 + 0j])
        for l in lam:
            poly = np.convolve(poly, np.array([1 - l, l]))  # (1-l) + l s
        pk = poly.real[:6]
        P0 = float(rig[rf]["det"].mid())
    kap = float(varkappa_arb(arb(fmpq(*rf))).mid())
    se0, sek = np.sqrt(p0 * (1 - p0) / nsamp), counts.std() / np.sqrt(nsamp)
    hist = [np.mean(counts == k) for k in range(5)]
    data[f"mc r={rv}"] = dict(n=nsamp, P0=p0, Ek=ek, hist=hist, det_pk=list(pk[:5]))
    R.check("prop:fredholm", f"Monte Carlo ({nsamp} samples of q, degree {D}) at r={rv}: P(k_r=0) and E k_r against "
            "det(I-X) and varkappa (within 4 standard errors)", abs(p0 - P0) < 4 * se0 + 1e-9 and abs(ek - kap) < 4 * sek, MC,
            f"P0: MC {p0:.4f} +- {se0:.4f} vs {P0:.4f}; E k: MC {ek:.4f} +- {sek:.4f} vs {kap:.4f}; "
            f"P(k=j), j<5: MC {[round(h, 4) for h in hist]} vs det {[round(float(p), 4) for p in pk[:5]]}")

R.finish(extra=data)
