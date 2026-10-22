"""Part 1: symbolic identities of main/main.tex.

Each identity is checked symbolically (SymPy) where that is feasible, for all parameters up to a stated size;
where it is not (power series, large n, transcendental parameters), it is checked at random points in
high-precision arithmetic (mpmath, 50 digits).  Inequalities with explicit constants are tested at random
points; such a test can refute a constant but does not prove it.
"""
import itertools
import random

import mpmath as mp
import sympy as sp

from common import Report, SYMBOLIC, EXACT, RANDOM, NUMERIC, INTERVAL, MC

mp.mp.dps = 50
R = Report("part1_symbolic")
rng = random.Random(20261022)
I = sp.I


def is_zero(e):
    e = sp.expand(e)
    if e == 0:
        return True
    return sp.simplify(e) == 0


def heat(expr, alpha, var):
    """exp(-alpha d^2/dvar^2) expr for a polynomial expr (finite sum)."""
    out, d, i = 0, sp.expand(expr), 0
    while d != 0:
        out += (-alpha) ** i / sp.factorial(i) * d
        d = sp.diff(d, var, 2)
        i += 1
    return sp.expand(out)


def He(k, x):
    return sp.hermite_prob(k, x) if k >= 0 else 0


def h(k, x):
    return He(k, x) / sp.sqrt(sp.factorial(k)) if k >= 0 else 0


def U(l, y):
    """Chebyshev U_l(y) with U_{-1}=0, U_{-2}=-1."""
    if l == -1:
        return 0
    if l == -2:
        return -1
    return sp.chebyshevu(l, y)


# ----------------------------------------------------------------------------------------------
# 1. The heat flow and the heat-Hermite identity
# ----------------------------------------------------------------------------------------------
n, t, z, x = sp.symbols("n t z x", positive=True)
K = 14
ok = all(is_zero(heat(x ** k, sp.Rational(1, 2), x) - He(k, x)) for k in range(K + 1))
R.check("lem:heat(a)", f"He_k(x) = exp(-d^2/2) x^k for k<={K}", ok, SYMBOLIC)
ok = all(is_zero(heat((sp.sqrt(n) * z) ** k, t / (2 * n), z) / sp.sqrt(sp.factorial(k))
                 - t ** sp.Rational(k, 2) * h(k, sp.sqrt(n / t) * z)) for k in range(K + 1))
R.check("lem:heat(a)", f"exp(-t/(2n) d_z^2)(sqrt(n)z)^k/sqrt(k!) = t^(k/2) h_k(sqrt(n/t) z), symbolic n,t,z, k<={K}",
        ok, SYMBOLIC)

ok = True
for N in range(1, 10):
    eta = sp.symbols(f"eta0:{N + 1}")
    P = sum(eta[N - k] * (sp.sqrt(N) * z) ** k / sp.sqrt(sp.factorial(k)) for k in range(N + 1))
    lhs = heat(P, t / (2 * N), z)
    rhs = t ** sp.Rational(N, 2) * sum(eta[j] * t ** sp.Rational(-j, 2) * h(N - j, sp.sqrt(N / t) * z)
                                      for j in range(N + 1))
    ok &= is_zero(lhs - rhs)
R.check("eq:heat-hermite", "P_{n,t}(z) = t^(n/2) sum_j eta_j t^(-j/2) h_{n-j}(sqrt(n/t) z), symbolic eta,t,z, n<=9",
        ok, SYMBOLIC)

# eq:reduction: leading zeros of q removed
ok = True
for N, m0 in [(5, 1), (6, 2), (7, 3)]:
    eta = sp.symbols(f"e0:{N + 1}")
    eta = [0] * m0 + list(eta[m0:])
    P = sum(eta[N - k] * (sp.sqrt(N) * z) ** k / sp.sqrt(sp.factorial(k)) for k in range(N + 1))
    Np = N - m0
    Pt = [eta[m0 + j] for j in range(Np + 1)]
    Ptil = sum(Pt[Np - k] * (sp.sqrt(Np) * z) ** k / sp.sqrt(sp.factorial(k)) for k in range(Np + 1))
    lhs = heat(P, t / (2 * N), z)
    rhs = heat(Ptil, t / (2 * Np), z).subs(z, sp.sqrt(sp.Rational(N, Np)) * z)
    ok &= is_zero(lhs - rhs)
R.check("eq:reduction", "P_{n,t}(z) = P~_{n',t}(sqrt(n/n') z) when eta_0..eta_{m0-1}=0 (n,m0)=(5,1),(6,2),(7,3)",
        ok, SYMBOLIC)


# lem:heat(b): exp(-alpha d^2) does not increase the number of nonreal zeros (randomised, exact counts)
def nonreal_count(p):
    P = sp.Poly(p, x)
    return P.degree() - len(P.real_roots())


ok, trials, worst = True, 300, 0
for _ in range(trials):
    d = rng.randint(2, 8)
    coeffs = [rng.randint(-6, 6) for _ in range(d)] + [rng.choice([-3, -2, -1, 1, 2, 3])]
    p = sum(c * x ** i for i, c in enumerate(coeffs))
    alpha = sp.Rational(rng.randint(1, 40), rng.randint(1, 20))
    q = heat(p, alpha, x)
    Zp, Zq = nonreal_count(p), nonreal_count(q)
    same_lead = sp.Poly(q, x).LC() == sp.Poly(p, x).LC() and sp.Poly(q, x).degree() == d
    ok &= (Zq <= Zp) and same_lead
R.check("lem:heat(b)", f"Z(exp(-a d^2)p) <= Z(p), degree and leading coefficient kept: {trials} random integer "
        "polynomials of degree 2..8, random rational a (exact Sturm counts)", ok, EXACT)

# ----------------------------------------------------------------------------------------------
# 2. Christoffel-Darboux
# ----------------------------------------------------------------------------------------------
ok = True
for N in range(1, 13):
    X = sp.sqrt(N) * z
    lhs = sum(h(k, X) ** 2 for k in range(N))
    rhs = N * (h(N, X) ** 2 - z * h(N, X) * h(N - 1, X) + h(N - 1, X) ** 2)
    ok &= is_zero(lhs - rhs)
R.check("lem:cd", "sum_{k<n} h_k(x)^2 = n(h_n^2 - z h_n h_{n-1} + h_{n-1}^2), x = sqrt(n) z, n<=12", ok, SYMBOLIC)

# ----------------------------------------------------------------------------------------------
# 3. Continuants, defect expansion
# ----------------------------------------------------------------------------------------------
y = sp.symbols("y")


def Pi(m, bs, var=y):
    """characteristic polynomial of the m x m tridiagonal matrix with zero diagonal, off-diagonal squares bs."""
    a, b = sp.Integer(1), var  # Pi_0, Pi_1
    if m == 0:
        return a
    for k in range(2, m + 1):
        a, b = b, sp.expand(var * b - bs[k - 2] * a)
    return b


JMAX = 12
s = sp.symbols(f"s0:{JMAX + 1}", positive=True)  # s_r = sqrt(beta_r), s_0 = 1


def AB_free(jmax):
    """A_{n,j}, B_{n,j} from the backward recurrence with sqrt(beta_r) replaced by free symbols s_r (s_0=1)."""
    sv = [sp.Integer(1)] + list(s[1:])
    A, B = [sp.Integer(1), sp.Integer(0)], [sp.Integer(0), sp.Integer(1)]
    for j in range(2, jmax + 1):
        A.append(sp.together((y * A[j - 1] - sv[j - 2] * A[j - 2]) / sv[j - 1]))
        B.append(sp.together((y * B[j - 1] - sv[j - 2] * B[j - 2]) / sv[j - 1]))
    return A, B


A_free, B_free = AB_free(JMAX)
ok = True
for j in range(2, JMAX + 1):
    d = sp.prod(s[1:j])
    beta = [s[r] ** 2 for r in range(1, JMAX)]
    ok &= is_zero(sp.cancel(B_free[j] * d) - Pi(j - 1, beta[0:j - 2]))
    ok &= is_zero(sp.cancel(A_free[j] * d) + Pi(j - 2, beta[1:j - 2]))
R.check("lem:continuant", f"B_(n,j) = Pi_(j-1)(beta_1..beta_(j-2))/d_(n,j), A_(n,j) = -Pi_(j-2)(beta_2..beta_(j-2))/d_(n,j), "
        f"j<={JMAX}, with sqrt(beta_r) free symbols (stronger than the stated identity)", ok, SYMBOLIC)

ok = all(is_zero(Pi(m, [1] * m) - U(m, y / 2)) for m in range(0, 13))
R.check("eq:fixed-j", "Pi_m(y;1,...,1) = U_m(y/2), m<=12 (limits A_j=-U_{j-2}, B_j=U_{j-1})", ok, SYMBOLIC)


def hvals_mp(N, X):
    hs = [mp.mpf(1), X]
    for k in range(1, N):
        hs.append((X * hs[k] - mp.sqrt(k) * hs[k - 1]) / mp.sqrt(k + 1))
    return hs


def AB_mp(N, Y, jmax):
    A, B = [mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]
    for j in range(2, jmax + 1):
        c1, c2 = mp.sqrt(1 - mp.mpf(j - 2) / N), mp.sqrt(1 - mp.mpf(j - 1) / N)
        A.append((Y * A[j - 1] - c1 * A[j - 2]) / c2)
        B.append((Y * B[j - 1] - c1 * B[j - 2]) / c2)
    return A, B


err = 0
for N in [7, 18, 30]:
    for _ in range(3):
        Y = mp.mpc(rng.uniform(-3, 3), rng.uniform(-1, 1))
        hs = hvals_mp(N, mp.sqrt(N) * Y)
        A, B = AB_mp(N, Y, N)
        for j in range(N + 1):
            err = max(err, abs(hs[N - j] - A[j] * hs[N] - B[j] * hs[N - 1]) / (abs(A[j] * hs[N]) + abs(B[j] * hs[N - 1])))
R.check("eq:AB", "h_{n-j}(sqrt(n)z) = A_{n,j} h_n + B_{n,j} h_{n-1} for all 0<=j<=n, n in {7,18,30}, random complex z "
        "(error relative to |A h_n|+|B h_{n-1}|)", err < mp.mpf(10) ** -45, RANDOM, f"max rel. error {mp.nstr(err, 3)}")


def matchings(m):
    edges = list(range(1, m))  # edge e joins vertices e and e+1 (vertices 1..m)
    out = []
    for r in range(0, m // 2 + 1):
        for S in itertools.combinations(edges, r):
            if all(S[i + 1] - S[i] >= 2 for i in range(len(S) - 1)):
                out.append(S)
    return out


ok = True
for m in range(1, 10):
    dl = sp.symbols(f"d1:{m}") if m > 1 else ()
    lhs = Pi(m, [1 - dl[e] for e in range(m - 1)])
    rhs = 0
    for S in matchings(m):
        removed = set()
        for e in S:
            removed |= {e, e + 1}
        comps, cur = [], 0
        for v in range(1, m + 1):
            if v in removed:
                if cur:
                    comps.append(cur)
                cur = 0
            else:
                cur += 1
        if cur:
            comps.append(cur)
        rhs += sp.prod([dl[e - 1] for e in S]) * sp.prod([U(c, y / 2) for c in comps])
    ok &= is_zero(lhs - rhs)
R.check("lem:defect", "matching expansion of Pi_m(y;1-delta_1..1-delta_{m-1}), m<=9, symbolic delta,y", ok, SYMBOLIC)

# ----------------------------------------------------------------------------------------------
# 4. Second-order continuant expansion and the envelope constant F_1^*
# ----------------------------------------------------------------------------------------------
eps = sp.symbols("epsilon", positive=True)


def A1(j, yy):
    if j <= 1:
        return sp.Integer(0)
    return (-sum((l + 1) * U(l - 1, yy / 2) * U(j - 3 - l, yy / 2) for l in range(1, j - 2))
            - sp.Rational(j * (j - 1), 4) * U(j - 2, yy / 2))


def B1(j, yy):
    if j <= 1:
        return sp.Integer(0)
    return (sum(l * U(l - 1, yy / 2) * U(j - 2 - l, yy / 2) for l in range(1, j - 1))
            + sp.Rational(j * (j - 1), 4) * U(j - 1, yy / 2))


ok = True
J2 = 9
subs = {s[r]: sp.sqrt(1 - r * eps) for r in range(1, JMAX + 1)}
for j in range(0, J2 + 1):
    for X_free, Xj, X1 in [(A_free, -U(j - 2, y / 2), A1(j, y)), (B_free, U(j - 1, y / 2), B1(j, y))]:
        ser = sp.series(X_free[j].subs(subs), eps, 0, 2).removeO()
        ok &= is_zero(ser - Xj - eps * X1)
R.check("lem:secondorder(b)", f"A_(n,j) = A_j + A^(1)_j/n + O(n^-2), same for B, j<={J2} (series in 1/n)", ok, SYMBOLIC)

w = sp.symbols("w")
yw = w + 1 / w
ok = True
Yseq = [sp.expand(A1(j, yw) + w * B1(j, yw)) for j in range(0, 13)]
ok &= Yseq[0] == 0 and Yseq[1] == 0
for j in range(2, 13):
    ok &= is_zero(Yseq[j] - (yw * Yseq[j - 1] - Yseq[j - 2] + sp.Rational(1, 2) * ((j - 1) * w ** j + (j - 2) * w ** (j - 2))))
R.check("lem:secondorder(d)", "Y_j = y Y_{j-1} - Y_{j-2} + ((j-1)w^j+(j-2)w^{j-2})/2, Y_0=Y_1=0, j<=12, symbolic w",
        ok, SYMBOLIC)

jj = sp.symbols("j")
al2 = (w ** 2 + 1) / (2 * (w ** 2 - 1))
al1 = (-sp.Rational(1, 2) * (w ** 2 + 2) - al2) / (w ** 2 - 1)
Pj = lambda k: al2 * k * (k + 1) / 2 + al1 * k
ok = is_zero(sp.simplify(w ** jj * Pj(jj) - yw * w ** (jj - 1) * Pj(jj - 1) + w ** (jj - 2) * Pj(jj - 2)
                         - sp.Rational(1, 2) * ((jj - 1) * w ** jj + (jj - 2) * w ** (jj - 2))))
k1, k2 = sp.symbols("k1 k2")
sol = sp.solve([Pj(0) + k1 + k2, w * (Pj(1) + k1) + k2 / w], [k1, k2])
ok &= all(is_zero(sp.together(w ** j * (Pj(j) + sol[k1]) + sol[k2] * w ** (-j) - Yseq[j])) for j in range(13))
R.check("lem:envelope(b)", "particular solution P(j)=alpha_2 j(j+1)/2+alpha_1 j of the Y-recurrence, and "
        "Y_j = w^j(P(j)+k1)+k2 w^-j with Y_0=Y_1=0 (j<=12)", ok, SYMBOLIC)
th = sp.symbols("theta", positive=True)
ew = sp.exp(I * th)
ok = is_zero(sp.simplify((al2.subs(w, ew) + I / 2 * sp.cot(th)).rewrite(sp.exp)))
ok &= is_zero(sp.simplify(((al1 + al2).subs(w, ew) - 1 / (4 * sp.sin(th) ** 2)).rewrite(sp.exp)))
R.check("lem:envelope(b)", "at w=e^{i theta}: alpha_2 = -(i/2)cot theta, alpha_1+alpha_2 = 1/(4 sin^2 theta), "
        "hence F_1^* = gamma(1/(4w^2) - (i/4)cot(theta) rho)", ok, SYMBOLIC)


# numerical check of F_1(T*,y*) = sum_j a_j Y_j against F_1^* for a concrete power series
def chebU_mp(L, Yv):
    us = [mp.mpf(1), 2 * (Yv / 2)]
    for l in range(1, L):
        us.append(Yv * us[l] - us[l - 1])
    return us


def Ufun(us, l):
    return 0 if l == -1 else (-1 if l == -2 else us[l])


qco = None  # q(w) = (1+2w)(1+4w/3+20w^2/9)/(1-5w/6), the symbol of Part 3b
def q_coeffs_rational(Jn):
    num = [mp.mpf(1), mp.mpf(4) / 3 + 2, mp.mpf(20) / 9 + mp.mpf(8) / 3, mp.mpf(40) / 9]
    out = []
    for j in range(Jn):
        out.append(sum(num[i] * (mp.mpf(5) / 6) ** (j - i) for i in range(min(j, 3) + 1)))
    return out


uq = mp.mpc(mp.mpf(-3) / 10, mp.mpf(6) / 10)  # zero of 1+4w/3+20w^2/9
qf = lambda ww: (1 + 2 * ww) * (1 + mp.mpf(4) / 3 * ww + mp.mpf(20) / 9 * ww ** 2) / (1 - mp.mpf(5) / 6 * ww)
dq = lambda ww, k=1: mp.diff(qf, ww, k)
Ts = 1 / abs(uq) ** 2
ths = mp.arg(uq)
wst = mp.exp(1j * ths)
yst = 2 * mp.cos(ths)
gam = uq * dq(uq)
rho = uq * dq(uq, 2) / dq(uq)
Jn = 400
eta = q_coeffs_rational(Jn)
aT = [eta[j] * Ts ** (-mp.mpf(j) / 2) for j in range(Jn)]
us = chebU_mp(Jn + 2, yst)
F1 = 0
for j in range(2, Jn):
    b1 = sum(l * Ufun(us, l - 1) * Ufun(us, j - 2 - l) for l in range(1, j - 1)) + mp.mpf(j * (j - 1)) / 4 * Ufun(us, j - 1)
    a1 = -sum((l + 1) * Ufun(us, l - 1) * Ufun(us, j - 3 - l) for l in range(1, j - 2)) - mp.mpf(j * (j - 1)) / 4 * Ufun(us, j - 2)
    F1 += aT[j] * (a1 + wst * b1)
om = mp.sin(ths)
F1star = gam * (1 / (4 * om ** 2) - 1j / 4 * mp.cot(ths) * rho)
R.check("lem:envelope(b)", "F_1(T*,y*) = sum_j a_j(A1_j+w*B1_j) equals F_1^* for q=(1+2w)(1+4w/3+20w^2/9)/(1-5w/6)",
        abs(F1 - F1star) < mp.mpf(10) ** -30, RANDOM, f"|diff| = {mp.nstr(abs(F1 - F1star), 3)}")

# F_0 derivatives at the collision point (closed form of the proof of lem:envelope(b))
def F0(tt, yy):
    ww = (yy + 1j * mp.sqrt(4 - yy ** 2)) / 2
    a = lambda v: qf(v / mp.sqrt(tt))
    return a(ww) + (wst - ww) * (a(ww) - a(1 / ww)) / (ww - 1 / ww)


beta = -1j * gam / (2 * om)
e1 = abs(F0(Ts, yst))
e2 = abs(mp.diff(lambda tt: F0(tt, yst), Ts) + gam / (2 * Ts))
e3 = abs(mp.diff(lambda yy: F0(Ts, yy), yst) - beta)
R.check("lem:envelope(b)", "F_0(T*,y*)=0, d_t F_0 = -gamma/(2T*), d_y F_0 = beta (same q)",
        max(e1, e2, e3) < mp.mpf(10) ** -30, RANDOM, f"errors {mp.nstr(e1, 2)}, {mp.nstr(e2, 2)}, {mp.nstr(e3, 2)}")

# lem:envelope(c): the affine identity behind the window
g1, g2, r1, r2, Tsym, xi, tau = sp.symbols("g1 g2 r1 r2 T xi tau", real=True)
gm, rh = g1 + I * g2, r1 + I * r2
omg = sp.sin(th)
bet = -I * gm / (2 * omg)
Wsym = Tsym / (2 * omg ** 2)
tauc = Tsym / 2 * sp.cot(th) * r2
xir = -sp.Rational(1, 2) * sp.cos(th) * r1
zeta_ = omg * (xi - xir)
lam_ = (tau - tauc - Wsym) / (2 * Wsym)
F1s = gm * (1 / (4 * omg ** 2) - I / 4 * sp.cot(th) * rh)
ok = is_zero(sp.simplify(bet * (zeta_ - I * lam_) / omg - (bet * xi - gm * tau / (2 * Tsym) + F1s)))
R.check("lem:envelope(c)", "beta(zeta - i lambda)/omega = beta xi - gamma tau/(2T*) + F_1^* with the window variables "
        "(eq:window) and tau_c, W, xi_r of the Notation", ok, SYMBOLIC)

# ----------------------------------------------------------------------------------------------
# 5. The limiting matrix: det H_a = a(w) a(1/w)
# ----------------------------------------------------------------------------------------------
zz = w + 1 / w
Aj = lambda j: -U(j - 2, zz / 2)
Bj = lambda j: U(j - 1, zz / 2)
ok = all(is_zero(sp.together(Aj(j) + w * Bj(j) - w ** j)) for j in range(0, 14))
ok &= all(is_zero(Aj(j + 1) + Bj(j)) and is_zero(Bj(j + 1) - zz * Bj(j) - Aj(j)) for j in range(0, 13))
R.check("lem:limit", "A_j + w B_j = w^j, A_{j+1} = -B_j, B_{j+1} = z B_j + A_j (z=w+1/w), j<=13", ok, SYMBOLIC)
ok = True
for d in range(0, 9):
    a = sp.symbols(f"a0:{d + 1}")
    cA = sum(a[j] * Aj(j) for j in range(d + 1))
    cB = sum(a[j] * Bj(j) for j in range(d + 1))
    H = sp.zeros(2, 2)
    for j in range(d + 1):
        H += a[j] * sp.Matrix([[Aj(j), Bj(j)], [Aj(j + 1), Bj(j + 1)]])
    ok &= all(is_zero(sp.together(e)) for e in H - sp.Matrix([[cA, cB], [-cB, cA + zz * cB]]))
    aw = sum(a[j] * w ** j for j in range(d + 1))
    ok &= is_zero(sp.together(cA ** 2 + zz * cA * cB + cB ** 2 - aw * aw.subs(w, 1 / w)))
R.check("lem:limit", "H_a = [[A,B],[-B,A+zB]] and det H_a(z) = a(w)a(1/w) for polynomial a of degree <=8 "
        "with symbolic coefficients", ok, SYMBOLIC)

# power series symbol: S1 of Part 3a, a(w) = (1-2w+2w^2)(1+5w/3)/(1-10w/13), radius 1.3
def S1coef(Jn):
    num = [mp.mpf(1), mp.mpf(5) / 3 - 2, 2 - mp.mpf(10) / 3, mp.mpf(10) / 3]
    return [sum(num[i] * (mp.mpf(10) / 13) ** (j - i) for i in range(min(j, 3) + 1)) for j in range(Jn)]


aS1 = S1coef(800)
aS1f = lambda v: (1 - 2 * v + 2 * v ** 2) * (1 + mp.mpf(5) / 3 * v) / (1 - mp.mpf(10) / 13 * v)
err = 0
for _ in range(6):
    if _ < 3:
        ww = mp.exp(1j * rng.uniform(0, mp.pi))
    else:
        ww = mp.mpf(rng.choice([-1, 1]) * rng.uniform(0.85, 1.15))
    Z = ww + 1 / ww
    us = chebU_mp(805, Z)
    cA = aS1[0] - sum(aS1[j] * Ufun(us, j - 2) for j in range(2, 800))
    cB = sum(aS1[j] * Ufun(us, j - 1) for j in range(1, 800))
    err = max(err, abs(cA ** 2 + Z * cA * cB + cB ** 2 - aS1f(ww) * aS1f(1 / ww)))
R.check("lem:limit", "det H_a(w+1/w) = a(w)a(1/w) for the power series (1-2w+2w^2)(1+5w/3)/(1-10w/13), |w|=1 and "
        "real w with e^-k0<=|w|<=e^k0", err < mp.mpf(10) ** -30, RANDOM, f"max error {mp.nstr(err, 3)}")

# ----------------------------------------------------------------------------------------------
# 6. Colleague matrix
# ----------------------------------------------------------------------------------------------
ok = True
for N in range(1, 6):
    a = sp.symbols(f"a0:{N + 1}")
    Jm = sp.zeros(N, N)
    for k in range(1, N):
        Jm[k - 1, k] = Jm[k, k - 1] = sp.sqrt(k)
    bvec = sp.Matrix([a[N - i] for i in range(N)])  # (a_n, ..., a_1)
    e = sp.zeros(N, 1)
    e[N - 1] = 1
    C = Jm / sp.sqrt(N) - (e * bvec.T) / a[0]
    lhs = (z * sp.eye(N) - C).det()
    Q = sum(a[j] * h(N - j, sp.sqrt(N) * z) for j in range(N + 1))
    ok &= is_zero(sp.together(lhs - sp.sqrt(sp.factorial(N)) / (a[0] * sp.Integer(N) ** sp.Rational(N, 2)) * Q))
R.check("lem:colleague", "det(zI - C) = sqrt(n!)/(a_0 n^(n/2)) Q(sqrt(n) z), n<=5, symbolic a_j", ok, SYMBOLIC)
err = 0
for N in [12, 25]:
    a = [mp.mpf(rng.uniform(-1, 1)) for _ in range(N + 1)]
    C = mp.matrix(N, N)
    for k in range(1, N):
        C[k - 1, k] = C[k, k - 1] = mp.sqrt(k) / mp.sqrt(N)
    for i in range(N):
        C[N - 1, i] -= a[N - i] / a[0]
    for _ in range(2):
        Z = mp.mpc(rng.uniform(-2, 2), rng.uniform(-1, 1))
        hs = hvals_mp(N, mp.sqrt(N) * Z)
        Q = sum(a[j] * hs[N - j] for j in range(N + 1))
        lhs = mp.det(Z * mp.eye(N) - C)
        rhs = mp.sqrt(mp.factorial(N)) / (a[0] * mp.mpf(N) ** (mp.mpf(N) / 2)) * Q
        err = max(err, abs(lhs - rhs) / abs(rhs))
R.check("lem:colleague", "same identity for n=12,25 at random complex z, random a_j", err < mp.mpf(10) ** -35, RANDOM,
        f"max rel. error {mp.nstr(err, 3)}")

# ----------------------------------------------------------------------------------------------
# 7. The model f = zeta cos(zeta+c) + lambda sin(zeta+c) and Kepler's equation
# ----------------------------------------------------------------------------------------------
ze, la, c, sv = sp.symbols("zeta lambda c s", real=True)
f = ze * sp.cos(ze + c) + la * sp.sin(ze + c)
fz = sp.diff(f, ze)
S_ = ze + c
Jf = lambda u: u - sp.sin(u) * sp.cos(u)
Rm = sp.Matrix([[sp.cos(S_), sp.sin(S_)], [-sp.sin(S_), sp.cos(S_)]])
gam_s = sp.Matrix([sp.sin(S_) * sp.cos(S_), -sp.cos(S_) ** 2])
ok = all(is_zero(e) for e in sp.Matrix([f, fz]) - Rm * (sp.Matrix([ze, la]) - gam_s))
ok &= is_zero(sp.expand_trig(f ** 2 + fz ** 2 - ((Jf(S_) - c) ** 2 + (la + sp.cos(S_) ** 2) ** 2)))
R.check("lem:model(a)", "Phi_c = R(s)[(zeta,lambda) - gamma(s)] and |(f, d_zeta f)|^2 = (J(s)-c)^2 + (lambda+cos^2 s)^2, "
        "s = zeta + c", ok, SYMBOLIC)
sc = sp.symbols("s_c", real=True)
sub = {ze: sp.sin(sc) * sp.cos(sc), la: -sp.cos(sc) ** 2, c: Jf(sc)}
vals = [sp.simplify(e.subs(sub)) for e in [f, fz, sp.diff(f, ze, 2), sp.diff(f, ze, 3)]]
DPhi = sp.Matrix([[sp.diff(f, ze), sp.diff(f, la)], [sp.diff(fz, ze), sp.diff(fz, la)]]).subs(sub)
ok = vals[0] == 0 and vals[1] == 0 and sp.simplify(vals[2] + 2 * sp.sin(sc)) == 0 \
    and sp.simplify(vals[3] + 2 * sp.cos(sc)) == 0 and sp.simplify(DPhi.det() - 2 * sp.sin(sc) ** 2) == 0
R.check("lem:model(b)", "at (zeta_c,lambda_c)=(sin s_c cos s_c, -cos^2 s_c), c=J(s_c): f=f'=0, f''=-2 sin s_c, "
        "f'''=-2cos s_c, det D Phi = 2 sin^2 s_c", ok, SYMBOLIC)
ok = sp.simplify(sp.diff(Jf(sv), sv) - 2 * sp.sin(sv) ** 2) == 0 and sp.simplify(Jf(-sv) + Jf(sv)) == 0 \
    and sp.simplify(Jf(sv + sp.pi) - Jf(sv) - sp.pi) == 0 and Jf(sp.pi / 2) == sp.pi / 2
ok &= sp.simplify(1 + 2 * (-sp.cos(sv) ** 2) + sp.cos(2 * sv)) == 0
Sb = sp.symbols("S")
ok &= sp.simplify((2 * Jf(Sb / 2)) - (Sb - sp.sin(Sb))) == 0 and sp.simplify(2 + 2 * (-sp.cos(Sb / 2) ** 2) - (1 - sp.cos(Sb))) == 0
R.check("sec:firstorder", "J' = 2 sin^2, J odd, J(s+pi)=J(s)+pi, J(pi/2)=pi/2; 1+2 lambda_c = -cos 2s_c; "
        "(2c, 2+2 lambda_c) on the cycloid (S - sin S, 1 - cos S), S=2 s_c", ok, SYMBOLIC)
uu = sp.symbols("u", real=True)
Fh = sp.exp(uu) * f.subs(la, la + 2 * uu)
ok = sp.simplify(sp.diff(Fh, uu) + sp.diff(Fh, ze, 2)) == 0
R.check("lem:model", "e^u f_{lambda+2u,c} solves d_u F = -d_zeta^2 F (model is a backward-heat trajectory)", ok, SYMBOLIC)


def s_of_c(cv):
    cv = mp.mpf(cv)
    return mp.findroot(lambda u: u - mp.sin(u) * mp.cos(u) - cv, cv)


def Psi(cv):
    d = abs(cv - mp.pi * mp.nint(cv / mp.pi))
    return mp.findroot(lambda P: P - mp.sin(P) - 2 * d, (mp.mpf(0.0), mp.pi), solver="bisect") if d > 0 else mp.mpf(0)


err = 0
for _ in range(40):
    cv = mp.mpf(rng.uniform(-10, 10))
    scv = s_of_c(cv)
    ct = cv - mp.pi * mp.floor(cv / mp.pi + mp.mpf(1) / 2)  # in [-pi/2, pi/2)
    P = Psi(cv)
    err = max(err, abs(P - 2 * abs(s_of_c(ct))), abs(mp.cos(2 * scv) - mp.cos(P)),
              abs(mp.sin(2 * scv) - mp.sign(ct) * mp.sin(P)))
R.check("eq:kepler", "Psi(c) = 2|s_{c~}| solves Psi - sin Psi = 2 dist(c, pi Z); cos 2s_c = cos Psi, "
        "sin 2s_c = sgn(c~) sin Psi (40 random c)", err < mp.mpf(10) ** -25, RANDOM, f"max error {mp.nstr(err, 3)}")
ok = abs(-mp.cos(s_of_c(0)) ** 2 + 1) < 1e-40 and abs(-mp.cos(s_of_c(mp.pi / 2)) ** 2) < 1e-40
R.check("lem:model", "Dini endpoints: lambda_0 = -1, lambda_{pi/2} = 0", ok, RANDOM)


# lem:model(c): real-rootedness threshold, by counting zeros (argument principle) and real zeros (sign changes)
def model_counts(cv, lv, Nn=5):
    import numpy as np
    fn = lambda u: u * np.cos(u + cv) + lv * np.sin(u + cv)
    # real zeros in (-N pi - c, N pi - c): sign changes on a fine grid plus local extrema checks
    a, b = -Nn * np.pi - cv, Nn * np.pi - cv
    xs = np.linspace(a, b, 400001)
    vs = fn(xs)
    nreal = int(np.sum(np.sign(vs[1:]) != np.sign(vs[:-1])))
    # complex zeros in |zeta + c| < N pi by the argument principle on the circle
    th_ = np.linspace(0, 2 * np.pi, 200001)
    zc = -cv + Nn * np.pi * np.exp(1j * th_)
    fv = zc * np.cos(zc + cv) + lv * np.sin(zc + cv)
    ang = np.unwrap(np.angle(fv))
    ntot = int(round((ang[-1] - ang[0]) / (2 * np.pi)))
    return nreal, ntot


ok, rows = True, []
for _ in range(12):
    cv = rng.uniform(0.02, mp.pi - 0.02)
    lc = float(-mp.cos(s_of_c(cv)) ** 2)
    for dl in [-0.05, 0.05]:
        nr, nt = model_counts(float(cv), lc + dl)
        expect_real = 11 if dl > 0 else 9
        ok &= (nt == 11 and nr == expect_real)
        rows.append((round(cv, 3), dl, nr, nt))
R.check("lem:model(c)", "all zeros of f_{lambda,c} in |zeta+c|<5pi real iff lambda>=lambda_c: 11 zeros there "
        "(argument principle), 11 real for lambda=lambda_c+0.05, 9 real for lambda_c-0.05 (12 random c)", ok, NUMERIC)

# lem:model(d) constants, with interval arithmetic
from flint import arb, ctx
ctx.prec = 200
pi = arb.pi()
CJ = arb(2).sqrt() * (arb(3) * pi ** 2 / 2) ** (arb(1) / 6)
c_z1 = 2 * (pi ** 3 / 12) ** (arb(2) / 3)
c_z2 = (arb(3) / 2 * pi ** 2) ** (arb(1) / 3) * (pi ** 3 / 3) ** (arb(2) / 3)
c_l1 = (pi / 2) ** (arb(1) / 3) + 2 * (pi ** 3 / 12) ** (arb(1) / 3)
c_l2 = (1 + CJ) * (pi ** 3 / 3) ** (arb(1) / 3)
ok = (CJ < arb("2.22")) and (c_z1 < arb("11.7")) and (c_z2 < arb("11.7")) and (c_l1 < arb("7.1")) and (c_l2 < arb("7.1"))
R.check("lem:model(d)", "C_J = sqrt2 (3pi^2/2)^(1/6) < 2.22 and the four constants below 11.7 (zeta) and 7.1 (lambda)", ok,
        INTERVAL, f"C_J={CJ.mid().str(5)}, zeta: {c_z1.mid().str(5)}, {c_z2.mid().str(5)}; lambda: {c_l1.mid().str(5)}, "
                  f"{c_l2.mid().str(5)}")
ok = (2 - arb(26) / 36 > 0) and (arb(9) / 2 * 2 + 1 + 2 <= 5 * (2 + 1 + 2))
R.check("thm:firstorder(c)", "arithmetic of the double-zero bound: 5(2/s^2+1+2/s)+1 <= 26/s^2 for s<=1, 2-26/36>0",
        ok and all((5 * (2 / sg ** 2 + 1 + 2 / sg) + 1) <= 26 / sg ** 2 + 1e-12 for sg in [k / 1000 for k in range(1, 1001)]),
        NUMERIC)


def Phi_c(cv, zv, lv):
    S = zv + cv
    return mp.sqrt((zv * mp.cos(S) + lv * mp.sin(S)) ** 2 + ((1 + lv) * mp.cos(S) - zv * mp.sin(S)) ** 2)


worst = 0
for _ in range(4000):
    cv = mp.mpf(rng.uniform(-4, 4))
    scv = s_of_c(cv)
    sv_ = scv + mp.mpf(rng.uniform(-1, 1)) * mp.mpf(10) ** rng.uniform(-4, 0.5)
    zv = sv_ - cv
    lv = -mp.cos(sv_) ** 2 + mp.mpf(rng.uniform(-1, 1)) * mp.mpf(10) ** rng.uniform(-5, 0)
    e = Phi_c(cv, zv, lv)
    if e > 1:
        continue
    dc = abs(cv - mp.pi * mp.nint(cv / mp.pi))
    d = max(e, dc)
    sg = abs(mp.sin(scv))
    lc, zc = -mp.cos(scv) ** 2, mp.sin(scv) * mp.cos(scv)
    r = [abs(zv - zc) / min((mp.mpf(1.5) * mp.pi ** 2 * e) ** (mp.mpf(1) / 3), 11.7 * e * d ** (-mp.mpf(2) / 3)),
         abs(lv - lc) / min(e + 2.22 * e ** (mp.mpf(2) / 3), 7.1 * e * d ** (-mp.mpf(1) / 3)),
         (12 / mp.pi ** 3) * dc / sg ** 3 if sg > 0 else 0]
    if e < sg ** 3 / 4:
        r += [abs(zv - zc) / (2 * e / sg ** 2), abs(lv - lc) / (e * (1 + 2 / sg))]
    worst = max(worst, max(r))
R.check("lem:model(d)", "localisation bounds for |Phi_c| = eps <= 1 (4000 random points near the threshold curve)",
        worst <= 1, RANDOM, f"max ratio lhs/rhs = {mp.nstr(worst, 4)}")
worst = 0
for _ in range(2000):
    cv = rng.uniform(0, 4)
    side = rng.randint(0, 3)
    u_ = rng.uniform(-1, 1)
    zv, lv = [(u_ * 1, 0.5), (u_ * 1, -1.5), (1, -0.5 + u_), (-1, -0.5 + u_)][side]
    worst = max(worst, 1 / (2 * Phi_c(mp.mpf(cv), mp.mpf(zv), mp.mpf(lv))))
R.check("prop:threshold", "|Phi_c| >= 1/2 on the boundary of B_0=[-1,1]x[-3/2,1/2] (2000 random points)", worst <= 1,
        RANDOM, f"max of 1/(2|Phi|) = {mp.nstr(worst, 4)}")

# ----------------------------------------------------------------------------------------------
# 8. Local form of a collision
# ----------------------------------------------------------------------------------------------
v, aa, sg, x0 = sp.symbols("v a sigma x0")
ok = True
for k in range(0, 15):
    rhs = sum((-aa) ** i / sp.factorial(i) * sp.factorial(k) / sp.factorial(k - 2 * i) * w ** (k - 2 * i) for i in range(k // 2 + 1))
    ok &= is_zero(heat(w ** k, aa, w) - rhs)
    Hp = sum((-sp.Rational(1, 2)) ** i * sp.factorial(k) / (sp.factorial(i) * sp.factorial(k - 2 * i)) * v ** (k - 2 * i) for i in range(k // 2 + 1))
    Hm = sum((sp.Rational(1, 2)) ** i * sp.factorial(k) / (sp.factorial(i) * sp.factorial(k - 2 * i)) * v ** (k - 2 * i) for i in range(k // 2 + 1))
    ok &= is_zero(Hp - He(k, v)) and is_zero(Hm - sp.expand(I ** (-k) * He(k, I * v)))
    for epsv, Hk in [(1, Hp), (-1, Hm)]:
        lhs = heat((z - x0) ** k, epsv * sg ** 2 / 2, z).subs(z, x0 + sg * v)
        ok &= is_zero(lhs - sg ** k * Hk)
R.check("lem:localform", "exp(-a d^2)w^k formula; H^+_k = He_k, H^-_k(v) = i^-k He_k(iv); "
        "P_t(x0+sigma v) = sum b_k sigma^k H^eps_k(v) for t = t0 + eps n sigma^2 (k<=14)", ok, SYMBOLIC)
cs = sp.symbols("p0:7")
P0 = sum(cs[i] * z ** i for i in range(7))
Pt = heat(P0, t / (2 * n), z)
ok = is_zero(sp.diff(Pt, t) + sp.diff(Pt, z, 2) / (2 * n))
Pz = sp.symbols("Pz")
Jsym = sp.Matrix([[-sp.Symbol("P2") / (2 * n), Pz], [-sp.Symbol("P3") / (2 * n), sp.Symbol("P2")]])
ok &= is_zero(Jsym.det().subs(Pz, 0) + sp.Symbol("P2") ** 2 / (2 * n))
R.check("rem:localform", "d_t P_t = -P_t''/(2n) (degree-6 generic P) and Jacobian determinant -P''^2/(2n) at a double zero",
        ok, SYMBOLIC)

# ----------------------------------------------------------------------------------------------
# 9. varkappa(r), the intensity, the tail constant
# ----------------------------------------------------------------------------------------------
r_, rho_ = sp.symbols("r rho", positive=True)
kap = (1 / sp.pi) * ((1 + r_ ** 2) / (1 - r_ ** 2) * sp.atan(r_) - sp.atanh(r_))
ok = sp.simplify(sp.diff(kap, r_) - 4 / sp.pi * r_ * sp.atan(r_) / (1 - r_ ** 2) ** 2) == 0 and kap.subs(r_, 0) == 0
ok &= sp.simplify(4 / sp.pi * (sp.atan(r_) / (2 * (1 - r_ ** 2)) - (sp.atanh(r_) + sp.atan(r_)) / 4) - kap) == 0
R.check("lem:pairs", "varkappa(0)=0, varkappa' = (4/pi) r arctan r/(1-r^2)^2 and the integrated-by-parts form", ok, SYMBOLIC)
ser = sp.series(sp.pi * kap, r_, 0, 11).removeO()
ok = sp.expand(ser - (sp.Rational(4, 3) * (r_ ** 3 + r_ ** 5) + sp.Rational(152, 105) * (r_ ** 7 + r_ ** 9))) == 0
R.check("lem:pairs", "pi varkappa(r) = 4/3 (r^3+r^5) + 152/105 (r^7+r^9) + O(r^11); varkappa = 4/(3pi) r^3(1+r^2+O(r^4))",
        ok, SYMBOLIC, f"series: {ser}")
Xr, Yr = sp.symbols("X Y", real=True)
zc = Xr + I * Yr
Fk = lambda a_, b_: (a_ - b_) / (1 - a_ * b_) ** 2
rho1 = 1 / (sp.pi * I) / sp.Abs(1 - zc ** 2) * Fk(zc, sp.conjugate(zc))
rhoc = 2 * Yr / (sp.pi * sp.Abs(1 - zc ** 2) * (1 - Xr ** 2 - Yr ** 2) ** 2)
ok = sp.simplify(sp.expand(rho1 - rhoc)) == 0
uu2, th2 = sp.symbols("u theta2", real=True)
ok &= sp.simplify(sp.expand(sp.Abs(1 - rho_ ** 2 * sp.exp(2 * I * th2)) ** 2 - ((1 + rho_ ** 2) ** 2 - 4 * rho_ ** 2 * sp.cos(th2) ** 2)).rewrite(sp.cos)) == 0
R.check("lem:pairs", "first intensity = case n=1 of the Matsumoto-Shirai Pfaffian: rho_1 = 2 Im z/(pi|1-z^2|(1-|z|^2)^2); "
        "|1-rho^2 e^{2i theta}|^2 = (1+rho^2)^2 - 4 rho^2 cos^2 theta", ok, SYMBOLIC)
err = 0
for _ in range(5):
    rv = mp.mpf(rng.uniform(0.05, 0.99))
    lhs = mp.quad(lambda tt: mp.sin(tt) / abs(1 - rv ** 2 * mp.exp(2j * tt)), [0, mp.pi / 2, mp.pi])
    mid = 2 * mp.quad(lambda uv: 1 / mp.sqrt((1 + rv ** 2) ** 2 - 4 * rv ** 2 * uv ** 2), [0, 1])
    err = max(err, abs(lhs - 2 * mp.atan(rv) / rv), abs(mid - mp.asin(2 * rv / (1 + rv ** 2)) / rv))
    kk = 4 / mp.pi * mp.quad(lambda p: p * mp.atan(p) / (1 - p ** 2) ** 2, [0, rv])
    err = max(err, abs(kk - (1 / mp.pi) * ((1 + rv ** 2) / (1 - rv ** 2) * mp.atan(rv) - mp.atanh(rv))))
R.check("lem:pairs", "angular integral = (2/rho) arctan rho = (1/rho) arcsin(2rho/(1+rho^2)); radial integral = varkappa(r)",
        err < mp.mpf(10) ** -35, RANDOM, f"max error {mp.nstr(err, 3)}")
cc_, xx_ = sp.symbols("c x_", positive=True)
ok = sp.simplify(sp.integrate(4 * sp.sqrt(cc_ * xx_), (xx_, 0, cc_)) - sp.Rational(8, 3) * cc_ ** 2) == 0
ok &= sp.Rational(8, 3) * (1 / (2 * sp.pi)) * 1 == 4 / (3 * sp.pi)
bb = sp.symbols("b", positive=True)
ok &= sp.simplify(sp.Rational(8, 3) * (1 / (2 * bb)) ** 2 * bb ** 2 / 3 - sp.Rational(2, 9)) == 0
R.check("thm:tail", "area of D_0(c) is 8c^2/3; constant 8/3 f(0)^2 E eta^2 = 4/(3pi) (Gaussian), 2/9 (uniform)", ok, SYMBOLIC)
ok = (sp.Rational(15, 2) + sp.Rational(13, 2) == 14) and (2 * sp.Rational(15, 2) * 2 * sp.Rational(13, 2) <= 240) \
    and (4 * (5 * 4) <= 240)
# area of the rectangle (15|e|+10G)(13|e|+8G) <= 240 (|e|+G)^2
e2s, Gs = sp.symbols("e G", positive=True)
ok &= sp.Poly(sp.expand(240 * (e2s + Gs) ** 2 - (15 * e2s + 10 * Gs) * (13 * e2s + 8 * Gs)), e2s, Gs).coeffs() and \
    all(cf >= 0 for cf in sp.Poly(sp.expand(240 * (e2s + Gs) ** 2 - (15 * e2s + 10 * Gs) * (13 * e2s + 8 * Gs)), e2s, Gs).coeffs())
rr = sp.symbols("rr", positive=True)
ok &= all(sp.simplify(sp.Rational(1, 4) ** (j - 3) * sp.Rational(27, 8) * sp.Rational(3, 2) ** (j - 3)
                      - (sp.Rational(3, 2)) ** j * sp.Rational(1, 4) ** (j - 3)) == 0 for j in range(3, 10))
R.check("thm:tail", "box (eq:box) has area <= 240(|eta_2|+Gamma)^2; r^(j-2)|z|^j <= r 4^-(j-3) (27/8)(3/2)^(j-3) "
        "for r<=1/4, |z|<=3/2", bool(ok), SYMBOLIC)

# ----------------------------------------------------------------------------------------------
# 10. The joint law: density of tau_c - W cos Psi(U)
# ----------------------------------------------------------------------------------------------
Ps, Ws, ss = sp.symbols("Psi W s", positive=True)
ok = sp.simplify(sp.diff((Ps - sp.sin(Ps)) / 2, Ps) - (1 - sp.cos(Ps)) / 2) == 0
dens_s = ((1 - sp.cos(Ps)) / sp.pi / (Ws * sp.sin(Ps))).subs(Ps, sp.acos(-ss / Ws))
pW = 1 / (sp.pi * Ws) * sp.sqrt((Ws + ss) / (Ws - ss))
ok &= sp.simplify((dens_s ** 2 - pW ** 2)) == 0  # both positive on (-W,W)
sx = sp.symbols("sx", real=True)
pWx = 1 / (sp.pi * Ws) * sp.sqrt((Ws + sx) / (Ws - sx))
tot = sp.simplify(sp.integrate(((1 - sp.cos(Ps)) / sp.pi), (Ps, 0, sp.pi)))
mean = sp.simplify(sp.integrate(-Ws * sp.cos(Ps) * (1 - sp.cos(Ps)) / sp.pi, (Ps, 0, sp.pi)))
ok &= tot == 1 and sp.simplify(mean - Ws / 2) == 0
num_tot = mp.quad(lambda u: 1 / mp.pi * mp.sqrt((1 + u) / (1 - u)), [-1, 0, 1])
ok &= abs(num_tot - 1) < 1e-20
R.check("thm:jointlaw", "Psi(U), U uniform on [0,pi/2], has density (1-cos Psi)/pi; s=-W cos Psi has density "
        "(pi W)^-1 ((W+s)/(W-s))^(1/2), total mass 1, mean W/2", ok, SYMBOLIC)
import numpy as np
from scipy.optimize import brentq

Wv = 1.2590
Us = np.random.default_rng(7).uniform(0, np.pi / 2, 20000)
svals = np.array([-Wv * np.cos(brentq(lambda P: P - np.sin(P) - 2 * u, 0, np.pi)) for u in Us])
grid = np.linspace(-Wv, Wv, 401)[1:-1]
Fth = lambda s0: float(mp.quad(lambda u: 1 / (mp.pi * Wv) * mp.sqrt((Wv + u) / (Wv - u)), [-Wv, s0]))
ks = max(abs(np.mean(svals <= g) - Fth(g)) for g in grid[::20])
R.check("thm:jointlaw", "empirical CDF of -W cos Psi(U), 20000 samples, against the integrated density (KS distance)",
        ks < 0.02, MC, f"KS = {ks:.4f}")

# equidistribution of the first-order law along n for the benchmark (formula only)
th_b, a_b = 1.1, 2.0
nu = (th_b - np.sin(th_b) * np.cos(th_b)) / np.pi
rho_b = 1 - 1j / np.tan(th_b)
u_b = a_b ** -0.5 * np.exp(1j * th_b)
gam_b = u_b * (2j * np.sqrt(a_b) * np.sin(th_b))
W_b = a_b / (2 * np.sin(th_b) ** 2)
tc_b = -a_b / 2 / np.tan(th_b) ** 2
cinf = np.angle(gam_b) - np.pi / 4 - th_b / 2 - 0.25 * np.sin(2 * th_b) * rho_b.real
ns = np.arange(1, 200001)
cn = cinf - ns * np.pi * nu
dist = np.abs(cn - np.pi * np.round(cn / np.pi))
from scipy.special import ellipkinc  # noqa (scipy present)
psi_tab_u = np.linspace(0, np.pi / 2, 20001)
psi_tab = np.array([brentq(lambda P: P - np.sin(P) - 2 * u, 0, np.pi) if u > 0 else 0 for u in psi_tab_u])
Fvals = tc_b - W_b * np.cos(np.interp(dist, psi_tab_u, psi_tab))
ks = max(abs(np.mean(Fvals - tc_b <= g * W_b / Wv) - Fth(g)) for g in grid[::20])
R.check("thm:jointlaw(a)", "for the benchmark (a=2, theta*=1.1) the values tau_c - W cos Psi(c_n^nat), n<=2e5, are "
        "distributed with density p_W(tau - tau_c) (KS distance); mean", ks < 0.01, NUMERIC,
        f"KS = {ks:.5f}, mean = {Fvals.mean():.5f} vs tau_c+W/2 = {tc_b + W_b / 2:.5f}")

# ----------------------------------------------------------------------------------------------
# 11. Notation, phase and carrier identities of Section 10
# ----------------------------------------------------------------------------------------------
Rs = sp.symbols("R", positive=True)
us_ = Rs * sp.exp(I * th)
Tst = Rs ** -2
ok = sp.simplify(Tst / 2 * sp.cot(th) * r2 - sp.re(us_) * r2 / (2 * Rs ** 2 * sp.im(us_))) == 0
ok &= sp.simplify(Tst / (2 * sp.sin(th) ** 2) - 1 / (2 * sp.im(us_) ** 2)) == 0
ok &= sp.simplify(sp.sqrt(Tst) * (-sp.cos(th) * r1 / 2) + sp.re(us_) * r1 / (2 * Rs ** 2)) == 0
ok &= sp.simplify(sp.sqrt(Tst) / (2 * sp.sin(th)) - 1 / (2 * sp.im(us_))) == 0
ok &= sp.simplify(2 * sp.cos(th) / (2 * sp.sqrt(Tst)) - sp.re(us_)) == 0
nn, taun, xin = sp.symbols("n tau_n xi_n", positive=True)
Lam = Tst + taun / nn
lhs = nn * (sp.sqrt(Lam) * (2 * sp.cos(th) + xin / nn) - sp.sqrt(Tst) * 2 * sp.cos(th))
rhs = sp.sqrt(Lam) * xin + 2 * sp.cos(th) * taun / (sp.sqrt(Lam) + sp.sqrt(Tst))
ok &= sp.simplify(lhs - rhs) == 0
ok &= sp.simplify(2 * Tst / (2 * sp.sin(th) ** 2) - 1 / sp.im(us_) ** 2) == 0
R.check("sec:firstorder", "tau_c and W in both forms; sqrt(T*) xi_r, sqrt(T*)/(2 omega), y*/(2 sqrt T*) in the heat variable; "
        "exact n(z_n - z*) formula; interval length 2W = 1/(Im u*)^2", ok, SYMBOLIC)

mm, xs_ = sp.symbols("M x_s", positive=True)  # M = m + 1/2
phim = sp.acos(xs_ / (2 * sp.sqrt(mm)))
gfun = mm * Jf(phim)
ok = sp.simplify(sp.diff(gfun, mm) - phim) == 0
ok &= sp.simplify((nn / 2 + sp.Rational(1, 4)) * (sp.sin(2 * th) - 2 * th) + 3 * sp.pi / 4 - (3 * sp.pi / 4 - (nn + sp.Rational(1, 2)) * Jf(th))) == 0
Th = sp.symbols("Theta", real=True)
ok &= sp.simplify(sp.expand_trig(sp.cos(Th) * sp.cos(th) - sp.cos(Th - th) + sp.sin(Th) * sp.sin(th))) == 0
argg = sp.symbols("arg_gamma", real=True)
cnat_intro = argg - sp.pi / 4 - (nn + sp.Rational(1, 2)) * th + nn / 2 * sp.sin(2 * th) - sp.sin(2 * th) * r1 / 4
cinf_ = argg - sp.pi / 4 - th / 2 - sp.sin(2 * th) * r1 / 4
ok &= sp.simplify(cnat_intro - (cinf_ - nn * Jf(th))) == 0
argc = sp.pi / 4 - (nn + sp.Rational(1, 2)) * th + nn / 2 * sp.sin(2 * th)
ok &= sp.simplify(argc + argg - sp.pi / 2 - sp.sin(2 * th) * r1 / 4 - cnat_intro) == 0
R.check("prop:phase", "g'(m) = phi_m; Szego's phase = 3pi/4 - (n+1/2)J(phi); cos(T)cos(th)-cos(T-th) = -sin T sin th; "
        "c_n^nat of the introduction = c_inf - n pi nu = c_n with arg c_n(x*) replaced by its asymptotics", ok, SYMBOLIC)

hn, hn1 = sp.symbols("h_n h_n1", real=True)
car = hn + I * (hn * sp.cos(th) - hn1) / sp.sin(th)
ok = sp.simplify(sp.re(car) - hn) == 0 and sp.simplify(sp.expand(sp.re(sp.exp(I * th) * car)) - hn1) == 0
ok &= sp.simplify(sp.sin(th) ** 2 * (sp.re(car) ** 2 + sp.im(car) ** 2) - (hn ** 2 - 2 * sp.cos(th) * hn * hn1 + hn1 ** 2)) == 0
M2 = sp.Matrix([[1, -sp.cos(th)], [-sp.cos(th), 1]])
ok &= set(sp.simplify(e) for e in M2.eigenvals().keys()) == {1 - sp.cos(th), 1 + sp.cos(th)}
R.check("lem:carrier(a)", "h_n = Re c_n, h_{n-1} = Re(w* c_n), omega^2|c_n|^2 = h_n^2 - y* h_n h_{n-1} + h_{n-1}^2; "
        "form eigenvalues 1 +- cos theta*", ok, SYMBOLIC)
ok = True
xi_ = sp.symbols("xi", real=True)
for N in range(2, 9):
    un = sp.exp(-x ** 2 / 4) * h(N, x)
    un1 = sp.exp(-x ** 2 / 4) * h(N - 1, x)
    ok &= sp.simplify(sp.diff(un, x, 2) + (N + sp.Rational(1, 2) - x ** 2 / 4) * un) == 0
    ok &= sp.simplify(sp.diff(un, x) - (sp.sqrt(N) * un1 - x / 2 * un)) == 0
    X_ = 2 * sp.sqrt(N) * sp.cos(th) + xi_ / sp.sqrt(N)
    Uxi = un.subs(x, X_)
    Upr = sp.diff(Uxi, xi_)
    carX = (h(N, x) + I * (h(N, x) * sp.cos(th) - h(N - 1, x)) / sp.sin(th)).subs(x, X_)
    ok &= sp.simplify(sp.exp(-X_ ** 2 / 4) * carX - (Uxi - I * Upr / sp.sin(th) - I * xi_ * Uxi / (2 * N * sp.sin(th)))) == 0
    epsx = (sp.Rational(1, 2) - xi_ * sp.cos(th) - xi_ ** 2 / (4 * N)) / N
    ok &= sp.simplify(sp.diff(Uxi, xi_, 2) + (sp.sin(th) ** 2 + epsx) * Uxi) == 0
R.check("lem:carrier(b)", "u_n''+(n+1/2-x^2/4)u_n=0, u_n'=sqrt(n)u_{n-1}-(x/2)u_n, eq:carrier-Z, and "
        "U''+(omega^2+eps(xi))U=0 with the stated eps (n<=8)", ok, SYMBOLIC)

# Plancherel-Rotach: arg and modulus of c_n(x*) (prop:phase), n up to 12800.  (Plain mpmath: naive ball
# arithmetic is useless here, its radii grow like e^{0.5 n} in the bulk.)
mp.mp.dps = 30
rows, worst_arg, worst_mod = [], 0, 0
for thv in [0.7, 1.1, 2.3]:
    for N in [200, 800, 3200, 12800]:
        xs = 2 * mp.sqrt(N) * mp.cos(thv)
        a0_, a1_ = mp.mpf(1), xs
        sc = mp.mpf(1)
        for k in range(1, N):
            a0_, a1_ = a1_, (xs * a1_ - mp.sqrt(k) * a0_) / mp.sqrt(k + 1)
        e = mp.exp(-xs * xs / 4)
        hN, hN1 = float(a1_ * e), float(a0_ * e)
        cr = complex(hN, (hN * np.cos(thv) - hN1) / np.sin(thv))
        pred = np.pi / 4 - (N + 0.5) * thv + N / 2 * np.sin(2 * thv)
        da = np.angle(cr * np.exp(-1j * pred))
        dm = abs(cr) / ((2 / (np.pi * N)) ** 0.25 * np.sin(thv) ** -0.5) - 1
        worst_arg, worst_mod = max(worst_arg, N * abs(da)), max(worst_mod, N * abs(dm))
        rows.append((thv, N, da, dm))
mp.mp.dps = 50
R.check("prop:phase", "arg c_n(x*) - (pi/4 - (n+1/2)theta + (n/2) sin 2theta) and |c_n(x*)|e^{-x*^2/4}/((2/(pi n))^(1/4) "
        "sin^(-1/2)) - 1 are O(1/n): n x error bounded for theta in {0.7,1.1,2.3}, n=200..12800", worst_arg < 1 and worst_mod < 1,
        NUMERIC, f"max n*|arg err| = {worst_arg:.3f}, max n*|mod err| = {worst_mod:.3f}")

# Mehler's formula at complex arguments
err = 0
for _ in range(4):
    rv = mp.mpf(rng.uniform(0.1, 0.8))
    X1, Y1 = mp.mpc(rng.uniform(-2, 2), rng.uniform(-2, 2)), mp.mpc(rng.uniform(-2, 2), rng.uniform(-2, 2))
    hx, hy = hvals_mp(400, X1), hvals_mp(400, Y1)
    ssum = sum(hx[k] * hy[k] * rv ** k for k in range(401))
    rhs = (1 - rv ** 2) ** -0.5 * mp.exp((2 * rv * X1 * Y1 - rv ** 2 * (X1 ** 2 + Y1 ** 2)) / (2 * (1 - rv ** 2)))
    err = max(err, abs(ssum - rhs) / abs(rhs))
R.check("lem:mehler", "Mehler's formula at random complex x,y and r in (0.1,0.8), series to k=400", err < mp.mpf(10) ** -30,
        RANDOM, f"max rel. error {mp.nstr(err, 3)}")

# ----------------------------------------------------------------------------------------------
# 12. The benchmark of rem:benchmark
# ----------------------------------------------------------------------------------------------
ab = sp.symbols("a", positive=True)
ub = ab ** sp.Rational(-1, 2) * sp.exp(I * th)
bq = -2 * sp.sqrt(ab) * sp.cos(th)
qb = lambda W_: 1 + bq * W_ + ab * W_ ** 2
ok = sp.simplify(sp.expand(qb(ub)).rewrite(sp.cos)) == 0
qpr = bq + 2 * ab * ub
ok &= sp.simplify((qpr - 2 * I * sp.sqrt(ab) * sp.sin(th)).rewrite(sp.cos)) == 0
rhob = sp.simplify(ub * 2 * ab / qpr)
ok &= sp.simplify((rhob - (1 - I * sp.cot(th))).rewrite(sp.exp)) == 0
tcb = (ab / 2) * sp.cot(th) * sp.im(1 - I * sp.cot(th))
Wb = ab / (2 * sp.sin(th) ** 2)
ok &= sp.simplify(tcb + ab / 2 * sp.cot(th) ** 2) == 0 and sp.simplify(tcb + Wb - ab / 2) == 0
ok &= sp.simplify(tcb - Wb + ab / 2 * (1 + sp.cos(th) ** 2) / sp.sin(th) ** 2) == 0
R.check("rem:benchmark", "q=1+bw+aw^2: q(u*)=0, q'(u*)=2i sqrt(a) sin theta, rho* = 1-i cot theta, tau_c=-(a/2)cot^2, "
        "W=a/(2 sin^2), tau_c+W=a/2, tau_c-W=-(a/2)(1+cos^2)/sin^2", ok, SYMBOLIC)
ok = True
bs_, as_ = sp.symbols("b a", real=True)
for N in range(2, 10):
    Qn = h(N, x) + bs_ * t ** sp.Rational(-1, 2) * h(N - 1, x) + as_ / t * h(N - 2, x)
    chi = N - 1 - as_ / t * sp.sqrt(N * (N - 1))
    rhs = (x + bs_ * sp.sqrt(N / t)) * He(N - 1, x) - chi * He(N - 2, x)
    ok &= is_zero(sp.sqrt(sp.factorial(N)) * Qn - rhs)
    # determinant representation with J'
    chis = sp.symbols("chi", positive=True)
    Jp = sp.zeros(N, N)
    for k in range(1, N):
        Jp[k - 1, k] = Jp[k, k - 1] = sp.sqrt(k)
    Jp[N - 1, N - 1] = -bs_ * sp.sqrt(N / t)
    Jp[N - 2, N - 1] = Jp[N - 1, N - 2] = sp.sqrt(chis)
    ok &= is_zero((x * sp.eye(N) - Jp).det() - ((x + bs_ * sp.sqrt(N / t)) * He(N - 1, x) - chis * He(N - 2, x)))
R.check("rem:benchmark", "sqrt(n!)Q_n^(t) = (x+b sqrt(n/t))He_{n-1} - chi He_{n-2}, and = det(x - J') (n<=9)", ok, SYMBOLIC)
ok = all(abs(Nv * (2 * mp.sqrt(mp.mpf(Nv) / (Nv - 1)) - 2) - 1) < 2.0 / Nv for Nv in range(2, 3000, 7))
b207 = 207 * (2 * mp.sqrt(mp.mpf(207) / 206) - 2)
R.check("rem:benchmark", "n(a sqrt(n/(n-1)) - a) = a/2 + O(1/n) (a=2); value 1.0036 at n=207",
        ok and abs(b207 - mp.mpf("1.0036")) < 5e-5, NUMERIC, f"207(2sqrt(207/206)-2) = {mp.nstr(b207, 6)}")
nub = (1.1 - np.sin(1.1) * np.cos(1.1)) / np.pi
ok = abs(nub - 0.2215) < 1e-4 and abs(tc_b + 0.2590) < 1e-4 and abs(W_b - 1.2590) < 1e-4
R.check("rem:benchmark", "a=2, theta*=1.1: nu = 0.2215..., tau_c = -0.2590..., W = 1.2590...", ok, NUMERIC,
        f"nu={nub:.6f}, tau_c={tc_b:.6f}, W={W_b:.6f}")

# ----------------------------------------------------------------------------------------------
# 13. Exterior-lemma algebra and Remark rem:full
# ----------------------------------------------------------------------------------------------
xv, dd = sp.symbols("x_v d", positive=True)
iv = sp.Symbol("i", integer=True, nonnegative=True)
ser_x = sp.series(xv * (1 + xv) / (1 - xv) ** 3, xv, 0, 30).removeO()
ok = sp.expand(ser_x - sum((k + 1) ** 2 * xv ** (k + 1) for k in range(29))) == 0
x7 = mp.mpf(7) / 10
ok &= abs(mp.nsum(lambda i: (i + 1) ** 2 * x7 ** (i + 1), [0, mp.inf]) - x7 * (1 + x7) / (1 - x7) ** 3) < 1e-35
W1, W2 = sp.symbols("W1 W2", real=True)
wc = W1 + I * W2
ok &= sp.simplify(sp.Abs(wc + 1 / wc - 2) ** 2 - (sp.Abs(wc - 1) ** 2 / sp.Abs(wc)) ** 2) == 0
C0s, R0s, dls, ns_ = sp.symbols("C0 R0 delta n", positive=True)
next_ = 2 * sp.sqrt(2) * C0s / (sp.sqrt(1 - R0s ** -2) * dls * dd ** sp.Rational(7, 2))
ok &= sp.simplify(2 * C0s / (sp.sqrt(1 - R0s ** -2) * dd ** 2 * ns_) * sp.sqrt(2) / dd ** sp.Rational(3, 2) - next_ * dls / ns_) == 0
ok &= sp.simplify(sp.Rational(1, 1) * (1 - (1 - dd) ** 2) - dd * (2 - dd)) == 0
R.check("lem:ext-quant", "sum (i+1)^2 x^(i+1) = x(1+x)/(1-x)^3; |J(w)-2| = |w-1|^2/|w|; the onset algebra giving n_ext",
        bool(ok), SYMBOLIC)
ok = all(mp.sqrt(mp.mpf(Nv) ** Nv / mp.factorial(Nv)) >= mp.e ** (mp.mpf(Nv) / 2) / (mp.e ** 2 * Nv) ** 0.25 for Nv in range(1, 400))
ok &= all(mp.factorial(Nv) <= mp.e * mp.sqrt(Nv) * (Nv / mp.e) ** Nv for Nv in range(1, 400))
R.check("rem:full", "n! <= e sqrt(n)(n/e)^n and (n^n/n!)^(1/2) >= e^(n/2)/(e^2 n)^(1/4), n<400", ok, RANDOM)

# ----------------------------------------------------------------------------------------------
# 14. Explicit-constant inequalities of Part I, at random points (refutation tests)
# ----------------------------------------------------------------------------------------------
def AB_real(N, Y, jmax):
    return AB_mp(N, mp.mpf(Y), jmax)


worst = {"chebyshev": 0, "transfer": 0, "head-rate(a)": 0, "head-rate(b)": 0, "sec(a)": 0, "sec(b)": 0, "sec(c)": 0}
for _ in range(300):
    kap_ = mp.mpf(rng.uniform(0.05, 2.0))
    Lz = 2 * mp.cosh(kap_)
    Y = mp.mpf(rng.uniform(-1, 1)) * Lz
    N = rng.choice([20, 50, 120, 400])
    jmax = N // 2
    A, B = AB_real(N, Y, jmax + 2)
    us = chebU_mp(jmax + 3, Y)
    for l in range(0, jmax):
        worst["chebyshev"] = max(worst["chebyshev"], abs(us[l]) / (mp.sinh((l + 1) * kap_) / mp.sinh(kap_)),
                                 (mp.sinh((l + 1) * kap_) / mp.sinh(kap_)) / ((l + 1) * mp.e ** (l * kap_)))
    ck = mp.coth(kap_) + 0.5
    k1, c1 = min(kap_, 1), 15 / min(kap_, 1)
    for j in range(0, jmax + 1):
        worst["transfer"] = max(worst["transfer"], (abs(A[j]) + abs(B[j])) / (2 * (j + 1) * mp.e ** (j * kap_ + ck * j ** 2 / N)))
        if (j + 1) ** 2 <= k1 * N and j >= 2:
            Aj_, Bj_ = -Ufun(us, j - 2), Ufun(us, j - 1)
            worst["head-rate(a)"] = max(worst["head-rate(a)"], (abs(A[j] - Aj_) + abs(B[j] - Bj_)) / (c1 * j ** 3 * mp.e ** ((j - 1) * kap_) / N))
        if (j + 2) ** 2 <= k1 * N:
            om_ = mp.sqrt(1 - mp.mpf(j) / N)
            Aj1, Bj1 = -Ufun(us, j - 1), Ufun(us, j)
            lhs = abs(om_ * A[j + 1] - Aj1) + abs(om_ * B[j + 1] - Bj1)
            if lhs > 0:
                worst["head-rate(b)"] = max(worst["head-rate(b)"], lhs / ((c1 + 3) * (j + 1) ** 3 * mp.e ** (j * kap_) / N))
# second-order lemma at complex y in D_kappa
for _ in range(150):
    kap_ = mp.mpf(rng.uniform(0.05, 1.5))
    rad = mp.e ** (kap_ * rng.uniform(-1, 1))
    wv = rad * mp.exp(1j * rng.uniform(0, 2 * mp.pi))
    Y = wv + 1 / wv
    us = chebU_mp(60, Y)
    for l in range(0, 50):
        worst["sec(a)"] = max(worst["sec(a)"], abs(us[l]) / (mp.sinh((l + 1) * kap_) / mp.sinh(kap_)))
    N = rng.choice([16, 64, 400, 2500])
    jm = int(mp.sqrt(N))
    A, B = AB_mp(N, Y, jm)
    for j in range(2, min(jm, 40) + 1):
        b1 = sum(l * Ufun(us, l - 1) * Ufun(us, j - 2 - l) for l in range(1, j - 1)) + mp.mpf(j * (j - 1)) / 4 * Ufun(us, j - 1)
        a1 = -sum((l + 1) * Ufun(us, l - 1) * Ufun(us, j - 3 - l) for l in range(1, j - 2)) - mp.mpf(j * (j - 1)) / 4 * Ufun(us, j - 2)
        worst["sec(b)"] = max(worst["sec(b)"], (abs(a1) + abs(b1)) / (j ** 4 * mp.e ** (j * kap_)))
        lhs = abs(A[j] + Ufun(us, j - 2) - a1 / N) + abs(B[j] - Ufun(us, j - 1) - b1 / N)
        worst["sec(c)"] = max(worst["sec(c)"], lhs / (4 * mp.e ** mp.coth(kap_) * j ** 6 * mp.e ** (j * kap_) / N ** 2))
for key, lab in [("chebyshev", "lem:chebyshev"), ("transfer", "lem:transfer"), ("head-rate(a)", "lem:head-rate(a)"),
                 ("head-rate(b)", "lem:head-rate(b)"), ("sec(a)", "lem:secondorder(a)"), ("sec(b)", "lem:secondorder(b)"),
                 ("sec(c)", "lem:secondorder(c)")]:
    R.check(lab, "explicit bound at random points (ratio lhs/rhs must be <= 1)", worst[key] <= 1, RANDOM,
            f"max ratio = {mp.nstr(worst[key], 4)}")

R.finish()
