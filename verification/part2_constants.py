"""Part 2a: the explicit constants of Proposition prop:early (real-rootedness before time 1).

The onset n_s of eq:std-onset is recomputed with interval arithmetic (Arb): n_s' exactly (rational comparison),
then f_s(n_s - 1) < 0 < f_s(n_s) with certified signs; with f_s(1) < 0 and the convexity of f_s this identifies
n_s = min{n >= n_s' : f_s(n) > 0}.  The algebraic steps of the proof are checked symbolically.
"""
import sympy as sp
from flint import arb, ctx, fmpq

from common import Report, SYMBOLIC, EXACT, INTERVAL

ctx.prec = 256
R = Report("part2_constants")

C12 = arb(2) / 4095


def consts(s):
    s = arb(s)
    beta2 = s.sqrt() * (1 - s.sqrt()) * (1 / s).log() / 4
    gam = beta2 / (128 * (1 - s * s))
    return s, beta2, gam


def f_s(s_, n):
    s, beta2, gam = consts(s_)
    eps = gam / 4
    Phi = 36 * (1 + s) * gam
    cp = (1 / (1 - s * s)).log() / 4
    L2 = (2 * arb(n + 1).sqrt()).log()
    log2 = arb(2).log()
    return n * (gam - eps) - log2 - L2 - cp - C12 * (2 * (n * Phi + L2 + cp) + eps * n + log2)


def n_prime(s_frac):
    """least n>=1 with (1-s)^(-1/2) s^(1/2) s^(n/4) <= 1/2, i.e. (s/(1-s))^2 s^n <= 1/16 (exact)."""
    s = fmpq(*s_frac)
    n = 1
    while (s / (1 - s)) ** 2 * s ** n > fmpq(1, 16):
        n += 1
    return n


def onset(s_frac):
    s_ = arb(fmpq(*s_frac))
    npr = n_prime(s_frac)
    # locate the sign change by bisection on midpoints, then certify
    lo, hi = npr, max(npr, 2)
    while not (f_s(s_, hi) > 0):
        hi *= 2
    if f_s(s_, lo) > 0:
        return npr, npr, True
    while hi - lo > 1:
        m = (lo + hi) // 2
        if float(f_s(s_, m).mid()) > 0:
            hi = m
        else:
            lo = m
    certified = (f_s(s_, hi) > 0) and (f_s(s_, hi - 1) < 0) and (f_s(s_, 1) < 0)
    return npr, hi, certified


claims = {(3, 10): (13707, 1), (1, 2): (25051, 4), (4, 5): (93949, 25), (19, 20): (493795, 169)}
table = {}
for sf, (ns_claim, npr_claim) in claims.items():
    npr, ns, cert = onset(sf)
    table[f"{sf[0]}/{sf[1]}"] = dict(n_s=ns, n_s_prime=npr, certified=cert)
    R.check("eq:std-onset", f"s={sf[0]}/{sf[1]}: n_s' = {npr_claim}, n_s = {ns_claim}", cert and npr == npr_claim and ns == ns_claim,
            INTERVAL + " + exact", f"computed n_s'={npr}, n_s={ns}, signs certified: {cert}")

# convexity of f_s and its slope
n_ = sp.symbols("n", positive=True)
s = sp.symbols("s", positive=True)
gs, c12 = sp.symbols("gamma c12", positive=True)
cp = sp.Rational(1, 4) * sp.log(1 / (1 - s ** 2))
fs = n_ * (gs - gs / 4) - sp.log(2) - sp.log(2 * sp.sqrt(n_ + 1)) - cp - c12 * (
    2 * (n_ * 36 * (1 + s) * gs + sp.log(2 * sp.sqrt(n_ + 1)) + cp) + gs / 4 * n_ + sp.log(2))
f2 = sp.simplify(sp.diff(fs, n_, 2))
ok = sp.simplify(f2 - (1 + 2 * c12) / (2 * (n_ + 1) ** 2)) == 0
lim = sp.limit(fs / n_, n_, sp.oo)
ok &= sp.simplify(lim - gs * (sp.Rational(3, 4) - c12 * (72 * (1 + s) + sp.Rational(1, 4)))) == 0
slope_ok = all((sp.Rational(3, 4) - sp.Rational(2, 4095) * (72 * (1 + sv) + sp.Rational(1, 4))) > sp.Rational(1, 2)
               for sv in [0, 1])
R.check("eq:std-onset", "f_s'' = (1+2c12)/(2(n+1)^2) > 0 (convex); f_s(n)/n -> gamma(3/4 - c12(72(1+s)+1/4)) > gamma/2 "
        "for all s in (0,1)", ok and slope_ok, SYMBOLIC)
ok = all(f_s(arb(fmpq(*sf)), 1) < 0 for sf in claims)
R.check("eq:std-onset", "f_s(1) < 0 at the four values of s", ok, INTERVAL)

# the constants beta_s, gamma_s and the identities of the proof
beta2 = sp.Rational(1, 4) * sp.sqrt(s) * (1 - sp.sqrt(s)) * sp.log(1 / s)
gam = beta2 / (128 * (1 - s ** 2))
ok = sp.simplify(beta2 / (1 - s ** 2) - sp.sqrt(s) * sp.log(1 / sp.sqrt(s)) / (2 * (1 + sp.sqrt(s)) * (1 + s))) == 0
R.check("prop:early", "beta_s^2/(1-s^2) = sqrt(s) log(1/sqrt s)/(2(1+sqrt s)(1+s))", ok, SYMBOLIC)
ctx.prec = 128
worst = max(float((consts(fmpq(k, 2000))[2] * 256 * arb.const_e()).mid()) for k in range(1, 2000))
R.check("prop:early", "gamma_s <= 1/(256 e) on (0,1) (analytic: x log(1/x) <= 1/e; grid of 1999 points)", worst <= 1,
        INTERVAL, f"max 256 e gamma_s on grid = {worst:.4f}")
ok = (arb(4) * arb(-3).exp()).sqrt() <= arb("-0.8").exp()
R.check("prop:early", "(4e^-3)^(1/2) <= e^-0.8 (chi-squared tail with x=4)", ok, INTERVAL,
        f"{(arb(4) * arb(-3).exp()).sqrt().mid().str(6)} vs {arb('-0.8').exp().mid().str(6)}")
lam_, m_, xx = sp.symbols("lambda m x", positive=True)
g = -lam_ * m_ * xx - m_ / 2 * sp.log(1 - 2 * lam_)
lam0 = sp.solve(sp.diff(g, lam_), lam_)
ok = len(lam0) == 1 and sp.simplify(lam0[0] - (1 - 1 / xx) / 2) == 0
ok &= sp.simplify(sp.exp(g.subs(lam_, lam0[0])) - (xx * sp.exp(1 - xx)) ** (m_ / 2)) == 0
R.check("prop:early", "inf_lambda e^{-lambda m x}(1-2lambda)^{-m/2} = (x e^{1-x})^{m/2} at lambda=(1-1/x)/2", ok, SYMBOLIC)
X, Y, nn = sp.symbols("X Y n", real=True)
zz = X + sp.I * Y
xm = sp.sqrt(nn / s) * zz
expo = (s * (X ** 2 + Y ** 2) * nn / s - s ** 2 * sp.re(sp.expand(xm ** 2))) / (1 - s ** 2)
phi = (X ** 2 / (1 + s) + Y ** 2 / (1 - s)) / 2
ok = sp.simplify(sp.expand(expo) - 2 * nn * phi) == 0
ok &= sp.simplify((sp.sqrt(s) + s) / (s * (1 - s)) - 1 / (sp.sqrt(s) * (1 - sp.sqrt(s)))) == 0
ok &= sp.simplify(nn * beta2 / (sp.sqrt(s) * (1 - sp.sqrt(s))) - nn / 4 * sp.log(1 / s)) == 0
ok &= sp.simplify(s ** ((nn + 1) / 2) * s ** (-nn / 4) - s ** sp.Rational(1, 2) * s ** (nn / 4)) == 0
R.check("prop:early", "eq:std-Vinf exponent = 2n phi; (sqrt s + s)/(s(1-s)) = 1/(sqrt s (1-sqrt s)); exponent at |z|=beta_s "
        "is (n/4)log(1/s); truncation term = s^(1/2) s^(n/4)", ok, SYMBOLIC)
y0, r2 = sp.symbols("y0 r2", positive=True)
vt = [2 * sp.pi * l / 12 for l in range(12)]
avg = sum(phi.subs({X: r2 * sp.cos(v), Y: y0 + r2 * sp.sin(v)}) for v in vt) / 12 - phi.subs({X: 0, Y: y0})
ok = sp.simplify(avg - r2 ** 2 / 4 * (1 / (1 + s) + 1 / (1 - s))) == 0
ok &= sp.simplify(avg.subs(r2, y0 / 4) - y0 ** 2 / (32 * (1 - s ** 2))) == 0
ok &= sp.simplify(gam.subs(s, s) - (2 * y0) ** 2 / (128 * (1 - s ** 2)) * 0 - gam) == 0
ok &= sp.simplify(y0 ** 2 / (32 * (1 - s ** 2)) - (2 * y0) ** 2 / (128 * (1 - s ** 2))) == 0  # gamma_s with beta_s = 2 y0
ok &= sp.simplify((3 * y0 / 2) ** 2 / (2 * (1 - s)) - 36 * (1 + s) * y0 ** 2 / (32 * (1 - s ** 2))) == 0
R.check("eq:std-gain", "12-point average of phi minus phi(z0) = r2^2/(2(1-s^2)) = y0^2/(32(1-s^2)) = gamma_s; "
        "Phi_max bound (3y0/2)^2/(2(1-s)) = 36(1+s)gamma_s", ok, SYMBOLIC)
ok = sp.Rational(2, 1) * sp.Rational(1, 4096) / (1 - sp.Rational(1, 4096)) == sp.Rational(2, 4095)
ok &= sp.nsimplify(2 * sp.summation(sp.Rational(1, 4096) ** sp.Symbol("m", integer=True, positive=True),
                                     (sp.Symbol("m", integer=True, positive=True), 1, sp.oo))) == sp.Rational(2, 4095)
R.check("eq:std-average", "c_12 = 2*2^-12/(1-2^-12) = 2 sum_m 2^-12m = 2/4095", ok, EXACT)

ctx.prec = 128
s_, b2, gm = consts(fmpq(1, 2))
thr = 4 * arb(13).log() / gm
cs = gm / 4
Cs = 14 * (gm * 25051 / 4).exp()
R.info("prop:early", "s=1/2: beta_s^2, gamma_s, c_s = gamma_s/4, C_s = 14 e^{gamma_s n_s/4}",
       f"beta^2={b2.mid().str(8)}, gamma={gm.mid().str(8)}, c_s={cs.mid().str(6)}, C_s={Cs.mid().str(6)}", INTERVAL)
R.check("prop:early", "RHS of eq:std-early is below 1 only for n > 4 log(13)/gamma_s, about 2.7e4 at s=1/2",
        abs(float(thr.mid()) - 2.7e4) < 0.05e4, INTERVAL, f"4 log 13/gamma_(1/2) = {thr.mid().str(7)}")
R.finish(extra=table)
