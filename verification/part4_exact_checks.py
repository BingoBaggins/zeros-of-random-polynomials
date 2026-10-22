"""Part 4: the exact rational-arithmetic computations of Remark rem:bernoulli-general (Remark 11.3).

Ported from exact_checks.py (parts A and B; part C concerned a finite model that is not in the paper).  All points are
Gaussian rationals, so every quantity lies in Q or Q(i) and no floating point enters a decision.

With c_j = 1/(pi |1 - z_j^2|), eq:ms-kernel gives rho_n(z_1..z_n) = prod c_j * i^{-n} Pf[K^c(z_i,z_j)], and
i^{-n} Pf[K^c] is real and rational at Gaussian-rational points.

(A) Hermitian determinantal test: T^2 <= 4 a12 a13 a23 for such a process; the c_j cancel in T^2/(4 a12 a13 a23).
(B) Four-disk counterexample: P(tau) = sum_S (-tau)^|S| prod_{j in S} a_j rho_|S| with a_j = pi b_j |1-z_j^2|, so
    a_j c_j = b_j and P has rational coefficients; its exact discriminant and Sturm count decide its real roots.

Cross-checks not in the original program: Pf(B)^2 = det(B) for every Pfaffian used (determinant by SymPy over Q(i)),
antisymmetry of the 2n x 2n matrix, rho_2(z,z) = 0, and the printed digits of the paper.
"""
from fractions import Fraction as Fr
from itertools import combinations

import mpmath as mp
import sympy as sp

from common import Report, EXACT, NUMERIC

R = Report("part4_exact_checks")
LABEL = "rem:bernoulli-general"


class GQ:
    """Gaussian rational re + i im, exact."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=0):
        self.re = Fr(re); self.im = Fr(im)

    def __add__(s, o): o = o if isinstance(o, GQ) else GQ(o); return GQ(s.re + o.re, s.im + o.im)
    __radd__ = __add__

    def __sub__(s, o): o = o if isinstance(o, GQ) else GQ(o); return GQ(s.re - o.re, s.im - o.im)

    def __rsub__(s, o): return GQ(o) - s

    def __neg__(s): return GQ(-s.re, -s.im)

    def __mul__(s, o):
        o = o if isinstance(o, GQ) else GQ(o)
        return GQ(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    __rmul__ = __mul__

    def __eq__(s, o): o = o if isinstance(o, GQ) else GQ(o); return s.re == o.re and s.im == o.im

    def conj(s): return GQ(s.re, -s.im)

    def abs2(s): return s.re * s.re + s.im * s.im

    def __truediv__(s, o):
        o = o if isinstance(o, GQ) else GQ(o)
        d = o.abs2(); n = s * o.conj()
        return GQ(n.re / d, n.im / d)

    def sym(s): return sp.Rational(s.re.numerator, s.re.denominator) + sp.I * sp.Rational(s.im.numerator, s.im.denominator)

    def __repr__(s): return f"{s.re}+{s.im}i"


def F(x, y):
    d = 1 - x * y
    return (x - y) / (d * d)


def Kc(z, w):
    return [[F(z, w), F(z, w.conj())], [F(z.conj(), w), F(z.conj(), w.conj())]]


def big(zs):
    n = len(zs); B = [[GQ(0)] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            k = Kc(zs[i], zs[j])
            for a in range(2):
                for b in range(2):
                    B[2 * i + a][2 * j + b] = k[a][b]
    return B


def pf(B):
    """Pfaffian by expansion along the first row."""
    n = len(B)
    if n == 0:
        return GQ(1)
    tot = GQ(0)
    for j in range(1, n):
        if B[0][j] == 0:
            continue
        idx = [k for k in range(1, n) if k != j]
        term = B[0][j] * pf([[B[a][b] for b in idx] for a in idx])
        tot = tot + term if j % 2 == 1 else tot - term
    return tot


PF_CHECKS = []   # (antisymmetric?, Pf^2 == det?) for every Pfaffian evaluated


def rho_scaled(zs):
    """i^{-n} Pf[K^c]_S = rho_n / prod c_j; must be real."""
    B = big(zs)
    p = pf(B)
    antisym = all(B[i][j] == -B[j][i] for i in range(len(B)) for j in range(len(B)))
    det = sp.Matrix(len(B), len(B), lambda i, j: B[i][j].sym()).det(method="berkowitz")
    PF_CHECKS.append((antisym, sp.expand(p.sym() ** 2 - det) == 0))
    v = [GQ(1), GQ(0, -1), GQ(-1), GQ(0, 1)][len(zs) % 4] * p
    assert v.im == 0, "i^{-n} Pf is not real"
    return v.re


# ---------- (A) three-point test ----------
def three_point(zs):
    k = [rho_scaled([z]) for z in zs]
    A = {(i, j): k[i] * k[j] - rho_scaled([zs[i], zs[j]]) for i, j in combinations(range(3), 2)}
    P3 = rho_scaled(zs)
    T = P3 - k[0] * k[1] * k[2] + k[0] * A[(1, 2)] + k[1] * A[(0, 2)] + k[2] * A[(0, 1)]
    return T * T / (4 * A[(0, 1)] * A[(1, 2)] * A[(0, 2)]), A


zs3 = [GQ(Fr(-3, 20), Fr(1, 20)), GQ(Fr(9, 10), Fr(7, 50)), GQ(Fr(-93, 100), Fr(13, 100))]
ratio, A = three_point(zs3)
R.check(LABEL, "the three points lie in the upper half disk D^+", all(z.im > 0 and z.abs2() < 1 for z in zs3), EXACT)
R.check(LABEL, "pair terms a_ij = rho1 rho1 - rho2 are positive at the three points", all(v > 0 for v in A.values()),
        EXACT, {k: float(v) for k, v in A.items()})
R.check(LABEL, "T^2/(4 a12 a13 a23) > 1, so no Hermitian determinantal kernel", ratio > 1, EXACT,
        f"ratio = {ratio.numerator}/{ratio.denominator} = {float(ratio):.10f}")
R.check(LABEL, "the ratio is 3.9154... as printed", f"{float(ratio):.10f}".startswith("3.9154"), NUMERIC, f"{float(ratio):.10f}")


# ---------- (B) four-disk counterexample ----------
def local_poly(zs, bs):
    n = len(zs); coeff = [Fr(0)] * (n + 1)
    for m in range(n + 1):
        for S in combinations(range(n), m):
            w = Fr(1)
            for j in S:
                w *= bs[j]
            coeff[m] += (-1) ** m * w * rho_scaled([zs[j] for j in S])
    return coeff  # P(tau) = sum coeff[m] tau^m


z1 = GQ(Fr(-19, 25), Fr(1, 10)); z3 = GQ(Fr(27, 100), Fr(17, 20))
zs4 = [z1, GQ(-z1.re, z1.im), z3, GQ(-z3.re, z3.im)]          # z2 = -conj z1, z4 = -conj z3
bs = [Fr(42, 5), Fr(42, 5), Fr(9, 40), Fr(9, 40)]
c = local_poly(zs4, bs)
Psym = sum(sp.Rational(x.numerator, x.denominator) * sp.Symbol("tau") ** m for m, x in enumerate(c))
tau = sp.Symbol("tau")
P = sp.Poly(Psym, tau, domain="QQ")

R.check(LABEL, "the four points lie in D^+ and the configuration is symmetric under z -> -conj z",
        all(z.im > 0 and z.abs2() < 1 for z in zs4), EXACT)
R.check(LABEL, "rho_2(z,z) = 0 (correlation functions vanish on the diagonal)",
        all(rho_scaled([z, z]) == 0 for z in zs4), EXACT)
R.check(LABEL, "P has degree 4, constant term 1, rational coefficients", P.degree() == 4 and c[0] == 1, EXACT,
        "; ".join(f"tau^{m}: {x}" for m, x in enumerate(c)))
printed = [1, -38.031, 534.81, -3288.5, 7440.7]
digits = [None, 3, 2, 1, 1]
ok = all(round(float(x), d) == p for x, p, d in zip(c, printed, digits) if d is not None)
R.check(LABEL, "coefficients match the printed P = 1 - 38.031 tau + 534.81 tau^2 - 3288.5 tau^3 + 7440.7 tau^4", ok,
        NUMERIC, [round(float(x), 4) for x in c])
D = sp.discriminant(P)
R.check(LABEL, "discriminant of P is negative (exact)", D < 0, EXACT, f"{float(D):.6e}")
R.check(LABEL, "P is squarefree: gcd(P, P') = 1 (exact)", sp.gcd(P, P.diff(tau)).degree() == 0, EXACT)
nreal = P.count_roots()   # Sturm, exact over Q
R.check(LABEL, "P has exactly two real roots (Sturm), hence one pair of simple nonreal roots", nreal == 2, EXACT,
        f"real roots: {nreal}")
mp.mp.dps = 30
roots = mp.polyroots([mp.mpf(x.numerator) / x.denominator for x in reversed(c)], maxsteps=200, extraprec=200)
nonreal = sorted((r for r in roots if abs(mp.im(r)) > mp.mpf(10) ** -20), key=lambda r: mp.im(r))
ok = len(nonreal) == 2 and all(abs(mp.re(r) - mp.mpf("0.08879")) < 5e-6 and abs(abs(mp.im(r)) - mp.mpf("0.00176")) < 5e-6
                               for r in nonreal)
R.check(LABEL, "the nonreal roots are 0.08879 +- 0.00176 i as printed", ok, NUMERIC,
        ", ".join(mp.nstr(r, 8) for r in nonreal))

R.check(LABEL, f"every Pfaffian used ({len(PF_CHECKS)}) is of an antisymmetric matrix and satisfies Pf^2 = det",
        all(a and b for a, b in PF_CHECKS), EXACT)
R.finish()
