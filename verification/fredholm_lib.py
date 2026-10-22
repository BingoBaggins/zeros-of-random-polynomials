"""Fredholm determinant det(I - X(r)) = P(T_* <= r^-2) of Proposition prop:fredholm, three ways.

(1) rigorous_det(r): Arb enclosure.  The entries G_mn(r) are given by the absolutely convergent series
        G_mn = (2/pi) sum_s T_k(s) r^(p+2s)/(p+2s),   k = m-n (odd), p = m+n+2,
        T_k(s) = sum_{a+b=s} c_a c_b 2k/(k^2 - 4(a-b)^2),  c_a = binom(2a,a)/4^a,
    obtained from |1 - x e^{2i theta}|^-1 = sum_{a,b} c_a c_b x^(a+b) cos(2(a-b)theta) and
    int_0^pi sin(k theta) cos(2j theta) d theta = 2k/(k^2-4j^2) (k odd).  Since sum_{a+b=s} c_a c_b = 1 and
    |k^2 - 4j^2| >= 1, |T_k(s)| <= 2|k|, which bounds the series tail.  The matrix is truncated to even indices
    m, m' <= M; by the bound |X_mm'| <= (2m+1) r^(m+m'+1)/(1-r^2) of the proof of prop:fredholm, the discarded
    part E has trace norm <= tau_M (sum of these bounds outside the block), and
        |det(I-X)/det(I-X_M) - 1| <= exp(||(I-X_M)^-1|| tau_M) - 1
    (det(I-X) = det(I-X_M) det(I-(I-X_M)^-1 E), |det(I+F)-1| <= e^{||F||_1}-1).  ||(I-X_M)^-1|| is bounded by the
    Frobenius norm of the Arb inverse.
(2) nystrom_det(r, nr, nt): Bornemann's method for the operator L of eq:Lkernel, restricted to the even functions
    (same nonzero eigenvalues), i.e. the kernel 2L(z,w) on the half disk D_r^+, with a Gauss-Legendre product rule
    in polar coordinates.
(3) float_X(r, M): the matrix X_M in double precision, with S_k(rho) from the Fourier coefficients of
    |1 - x e^{i phi}|^-1 (FFT) and composite Gauss-Legendre quadrature in rho graded towards r.
"""
import numpy as np
from flint import arb, arb_mat, ctx, fmpq


# ----------------------------------------------------------------------------------------------
# (1) rigorous
# ----------------------------------------------------------------------------------------------
def outside_bound(r, M):
    """sum of (2m+1) r^(m+m'+1)/(1-r^2) over even m,m' >= 0 with max(m,m') > M (arb, exact closed forms)."""
    r2 = r * r
    Ainf = 4 * r2 / (1 - r2) ** 2 + 1 / (1 - r2)
    Binf = 1 / (1 - r2)
    AM = sum(((4 * k + 1) * r ** (2 * k) for k in range(M // 2 + 1)), arb(0))
    BM = sum((r ** (2 * k) for k in range(M // 2 + 1)), arb(0))
    return r / (1 - r2) * (Ainf * Binf - AM * BM)


def rigorous_det(r_frac, M, N, prec=192):
    """Arb enclosure of det(I - X(r)) and of tr X(r); r = fmpq(*r_frac)."""
    ctx.prec = prec
    r = arb(fmpq(*r_frac))
    kmax = M + 1
    c = [arb(1)]
    for a in range(1, N + 1):
        c.append(c[-1] * (2 * a - 1) / (2 * a))
    # T[k][s] for odd k in 1..kmax:  T_k(s) = 2 sum_a c_a c_{s-a} / (k - 2s + 4a)
    T = {}
    prods = [[c[a] * c[s - a] for a in range(s + 1)] for s in range(N + 1)]
    for k in range(1, kmax + 1, 2):
        row = []
        for s in range(N + 1):
            ps = prods[s]
            acc = arb(0)
            for a in range(s + 1):
                acc += ps[a] / (k - 2 * s + 4 * a)
            row.append(2 * acc)
        T[k] = row
    pw = [arb(1)]
    for e in range(1, 2 * M + 2 * N + 8):
        pw.append(pw[-1] * r)
    two_over_pi = 2 / arb.pi()
    r2 = r * r
    cache = {}

    def G(mp_, n_):
        key = (mp_, n_)
        if key in cache:
            return cache[key]
        k, p = mp_ - n_, mp_ + n_ + 2
        sgn = 1 if k > 0 else -1
        Tk = T[abs(k)]
        acc = arb(0)
        for s in range(N + 1):
            acc += Tk[s] * pw[p + 2 * s] / (p + 2 * s)
        tail = 2 * abs(k) * pw[p + 2 * N + 2] / ((p + 2 * N + 2) * (1 - r2))
        val = two_over_pi * (sgn * acc + arb(0, tail.upper()))
        cache[key] = val
        return val

    idx = list(range(0, M + 1, 2))
    Xm = arb_mat(len(idx), len(idx))
    for i, m in enumerate(idx):
        for jx, mq in enumerate(idx):
            v = -(m + 1) * G(mq, m + 1)
            if m > 0:
                v += m * G(mq, m - 1)
            Xm[i, jx] = v
    n = len(idx)
    IX = arb_mat(n, n)
    for i in range(n):
        for jx in range(n):
            IX[i, jx] = (1 if i == jx else 0) - Xm[i, jx]
    detM = IX.det()
    inv = IX.inv()
    fro = arb(0)
    for i in range(n):
        for jx in range(n):
            fro += inv[i, jx] ** 2
    opn = max(fro.sqrt().upper(), arb(1))
    tauM = outside_bound(r, M)
    eta = opn * tauM
    rel = (eta.exp() - 1).upper()
    det_encl = detM * (1 + arb(0, rel))
    trM = sum((Xm[i, i] for i in range(n)), arb(0))
    # diagonal tail for the trace: |X_mm| <= (2m+1) r^(2m+1)/(1-r^2), m > M even
    trtail = sum(((2 * m + 1) * r ** (2 * m + 1) / (1 - r2) for m in range(M + 2, 3 * M + 400, 2)), arb(0)) + r ** (3 * M + 400)
    tr_encl = trM + arb(0, (trtail * 2).upper())
    return dict(det=det_encl, detM=detM, tauM=tauM, invnorm=opn, rel=rel, trace=tr_encl, X=Xm)


def varkappa_arb(r):
    return ((1 + r * r) / (1 - r * r) * r.atan() - r.atanh()) / arb.pi()


# ----------------------------------------------------------------------------------------------
# (2) Bornemann / Nystrom on the half disk
# ----------------------------------------------------------------------------------------------
def nystrom_matrix(r, nr, nt):
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xt, wt = np.polynomial.legendre.leggauss(nt)
    rho = r * (xr + 1) / 2
    th = np.pi * (xt + 1) / 2
    RR, TT = np.meshgrid(rho, th, indexing="ij")
    W = (np.outer(wr * r / 2, wt * np.pi / 2) * RR).ravel()
    zz = (RR * np.exp(1j * TT)).ravel()
    zc = np.conj(zz)
    Z = zz[:, None]
    Wb = zc[None, :]
    K = (1 / (np.pi * 1j * np.abs(1 - zz ** 2)))[None, :] * ((Z - Wb) / (1 - Z * Wb) ** 2 - (Z + Wb) / (1 + Z * Wb) ** 2)
    sw = np.sqrt(W)
    return sw[:, None] * K * sw[None, :]


def nystrom_det(r, nr, nt):
    A = nystrom_matrix(r, nr, nt)
    sign, logdet = np.linalg.slogdet(np.eye(A.shape[0]) - A)
    return sign * np.exp(logdet), logdet, sign


# ----------------------------------------------------------------------------------------------
# (3) double-precision matrix X_M
# ----------------------------------------------------------------------------------------------
def graded_nodes(r, q=24, levels=None):
    if levels is None:
        levels = int(np.ceil(np.log2(r / max(1 - r, 1e-3)))) + 6
    edges = [0.0] + [r - r * 2.0 ** (-i) for i in range(1, levels + 1)] + [r]
    x, w = np.polynomial.legendre.leggauss(q)
    nodes, weights = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        nodes.append((b - a) / 2 * x + (a + b) / 2)
        weights.append((b - a) / 2 * w)
    return np.concatenate(nodes), np.concatenate(weights)


def float_X(r, M, q=24, L=None, levels=None):
    kmax = M + 1
    rho, wts = graded_nodes(r, q, levels)
    if L is None:
        L = int(2 ** np.ceil(np.log2(max(256, 60 / (1 - r * r)))))
    Jmax = L // 2
    js = np.arange(Jmax + 1)
    ks = np.arange(1, kmax + 1, 2)
    Kmat = 2 * ks[:, None] / (ks[:, None] ** 2 - 4.0 * js[None, :] ** 2)
    phi = 2 * np.pi * np.arange(L) / L
    mev = np.arange(0, M + 1, 2)
    nod = np.arange(1, M + 2, 2)
    kk = mev[:, None] - nod[None, :]
    kidx = (np.abs(kk) - 1) // 2
    sg = np.sign(kk)
    G = np.zeros((len(mev), len(nod)))
    for p_, w_ in zip(rho, wts):
        x = p_ * p_
        g = 1 / np.abs(1 - x * np.exp(1j * phi))
        F = np.fft.rfft(g).real / L
        gh = np.empty(Jmax + 1)
        gh[0] = F[0]
        gh[1:] = 2 * F[1:Jmax + 1]
        S = Kmat @ gh
        G += w_ * np.outer(p_ ** mev, p_ ** (nod + 1)) * (sg * S[kidx])
    G *= 2 / np.pi
    # X[m, m'] = m G[m', m-1] - (m+1) G[m', m+1]; G rows indexed by m' (even), columns by n (odd)
    n = len(mev)
    X = np.zeros((n, n))
    for i, m in enumerate(mev):
        col_p = (m + 1 - 1) // 2  # index of n = m+1
        X[i, :] = -(m + 1) * G[:, col_p]
        if m > 0:
            X[i, :] += m * G[:, (m - 1 - 1) // 2]
    return X


def logdet_cond(X):
    A = np.eye(X.shape[0]) - X
    sign, ld = np.linalg.slogdet(A)
    Ainv = np.linalg.inv(A)
    kappa = float(np.sum(np.abs(Ainv.T * X)))  # componentwise relative condition number of det(I-X)
    return sign, ld, kappa
