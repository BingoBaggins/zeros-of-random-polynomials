"""Certified counting of the real and nonreal zeros of Q(x) = sum_{j<=n} a_j h_{n-j}(x).

Certificate.  Let x_1 < ... < x_m be approximate real zeros and p_0 = -oo < p_1 < ... < p_{m-1} < p_m = +oo test
points between them (p_i exact dyadic numbers).  The sign of Q(p_i) is decided with Arb ball arithmetic (the ball must
exclude 0), the signs at -oo/+oo from the leading coefficient a_0/sqrt(n!).  Each sign change gives a real zero in a
distinct interval, so #real zeros >= #sign changes =: s.  For approximate nonreal zeros c (Im c > 0) a Krawczyk test
on a small box B around c with B cap R empty,
        K = c - Y Q(c) + (1 - Y Q'(B)) (B - c) contained in int(B),  0 not in Q'(B),  Y = 1/Q'(c),
proves a unique zero in B, and its conjugate gives a second nonreal zero (for real coefficients).  With pairwise
disjoint boxes this gives #nonreal >= 2k.  If s + 2k = n, both counts are exact (and the real zeros are simple).

Naive ball evaluation of the three-term recurrence over-estimates radii by roughly e^{0.55 n} in the bulk (the radii
obey the recurrence with absolute values), so the working precision is ~1.6 n + 200 bits and the Krawczyk boxes have
radius ~2^-(0.8n+60); approximate zeros are refined by Newton's method at that precision.
"""
import numpy as np
from flint import acb, arb, ctx


class HermiteCombination:
    """Q = sum_j a_j h_{n-j}; coeff_fn(prec) returns [a_0..a_n] as arb at the current precision."""

    def __init__(self, coeff_fn, n):
        self.coeff_fn = coeff_fn
        self.n = n
        self.prec = None

    def set_prec(self, prec):
        if self.prec == prec:
            return
        ctx.prec = prec
        self.prec = prec
        n = self.n
        self.a = self.coeff_fn(prec)
        self.isq = [1 / arb(k + 1).sqrt() for k in range(n + 1)]  # 1/sqrt(k+1)
        self.sq = [arb(k).sqrt() for k in range(n + 1)]
        # weights: Q = sum_k c0[k] h_k, Q' = sum_k c1[k] h_{k-1}, Q'' = sum_k c2[k] h_{k-2}
        self.c0 = [self.a[n - k] for k in range(n + 1)]
        self.c1 = [self.a[n - k] * self.sq[k] for k in range(n + 1)]
        self.c2 = [self.a[n - k] * self.sq[k] * self.sq[k - 1] if k >= 2 else arb(0) for k in range(n + 1)]
        self.floats = np.array([float(c.mid()) for c in self.a])

    def eval(self, x, deriv=1):
        """Q(x), Q'(x) (and Q''(x)) for arb/acb x by the forward recurrence."""
        n = self.n
        hm1, h0 = x * 0, x * 0 + 1  # h_{-1}, h_0
        Q = self.c0[0] * h0
        D1 = x * 0
        D2 = x * 0
        hkm2, hkm1 = hm1, h0  # h_{k-2}, h_{k-1} at k=1
        for k in range(1, n + 1):
            hk = (x * hkm1 - self.sq[k - 1] * hkm2) * self.isq[k - 1]
            Q += self.c0[k] * hk
            D1 += self.c1[k] * hkm1
            if deriv >= 2 and k >= 2:
                D2 += self.c2[k] * hkm2
            hkm2, hkm1 = hkm1, hk
        return (Q, D1, D2) if deriv >= 2 else (Q, D1)

    # ------------------------------------------------------------------ approximate zeros
    def approx_zeros(self):
        """eigenvalues of the colleague matrix (lem:colleague), in the variable x = sqrt(n) z."""
        n = self.n
        a = self.floats
        C = np.zeros((n, n))
        k = np.arange(1, n)
        C[k - 1, k] = C[k, k - 1] = np.sqrt(k)
        C[n - 1, :] -= np.sqrt(n) * a[n:0:-1] / a[0]
        return np.linalg.eigvals(C)

    def newton(self, x0, iters=60, tol_bits=None):
        """Newton refinement of a (complex) zero at the current precision, on midpoints."""
        x = acb(x0)
        tol = arb(2) ** (-(tol_bits or (self.prec - 40)))
        for _ in range(iters):
            Q, D1 = self.eval(acb(x.real.mid(), x.imag.mid()))
            step = Q / D1
            step = acb(step.real.mid(), step.imag.mid())
            x = x - step
            if abs(step).mid() < tol * (1 + abs(x).mid()):
                break
        return acb(x.real.mid(), x.imag.mid())

    def krawczyk(self, c, rad):
        """True if the box c +- rad (re and im) contains a unique zero and does not meet R."""
        c = acb(c.real.mid(), c.imag.mid())
        rad = arb(rad)
        if not (c.imag.mid() - rad > 0):
            return False
        B = acb(arb(c.real.mid(), rad), arb(c.imag.mid(), rad))
        Qc, D1c = self.eval(c)
        _, D1B = self.eval(B)
        if D1B.contains(0):
            return False
        Y = 1 / D1c
        Y = acb(Y.real.mid(), Y.imag.mid())
        K = c - Y * Qc + (1 - Y * D1B) * (B - c)
        inner = acb(arb(c.real.mid(), rad * arb("0.999")), arb(c.imag.mid(), rad * arb("0.999")))
        return bool(inner.contains(K))

    def sign_at(self, x):
        Q = self.eval(arb(x), deriv=1)[0]
        if Q > 0:
            return 1
        if Q < 0:
            return -1
        return 0


def classify(H, want_full=False, prec=None, max_tries=3, imag_tol=1e-9, susp=1e-3):
    """Certify the zero structure of H (HermiteCombination).

    Returns dict(status, n_real_lb, n_nonreal_lb, exact, ...):
      'all-real'  n certified sign changes (all zeros real and simple);
      'nonreal'   at least one certified nonreal zero (Krawczyk); 'exact' tells whether real + nonreal = n;
      'undecided' otherwise.
    If want_full, all nonreal candidates are certified (for exact counts); else one suffices.
    """
    n = H.n
    if prec is None:
        prec = int(1.6 * n) + 200
    res = dict(status="undecided")
    for attempt in range(max_tries):
        H.set_prec(prec)
        lead_sign = 1 if H.a[0] > 0 else (-1 if H.a[0] < 0 else 0)
        if lead_sign == 0:
            raise ValueError("a_0 not certified nonzero")
        z = H.approx_zeros()
        sc_ = 1 + np.abs(z)
        realc = sorted(z[np.abs(z.imag) <= imag_tol * sc_].real)
        cplx = [c for c in z if c.imag > imag_tol * (1 + abs(c))]
        # clusters: nearly real complex candidates, and close real pairs
        clusters = [c.real for c in cplx if c.imag < susp * (1 + abs(c))]
        cplx = [c for c in cplx if c.imag >= susp * (1 + abs(c))]
        if len(realc) > 2:
            g = np.diff(realc)
            med = np.median(g)
            close = set()
            for i in np.where(g < med / 20)[0]:
                if i in close or i + 1 in close:
                    continue
                close |= {i, i + 1}
                clusters.append((realc[i] + realc[i + 1]) / 2)
            realc = [x for i, x in enumerate(realc) if i not in close]
        undecided = False
        for x0 in clusters:
            xe = _critical_point(H, x0)
            Qe, D1e, D2e = H.eval(arb(xe), deriv=2)
            if (Qe > 0 and D2e < 0) or (Qe < 0 and D2e > 0):
                d = float(abs(2 * Qe / D2e).sqrt().mid())
                realc += [xe - d, xe + d]
            elif (Qe > 0 and D2e > 0) or (Qe < 0 and D2e < 0):
                cplx.insert(0, complex(xe, float(abs(2 * Qe / D2e).sqrt().mid())))
            else:
                undecided = True
        realc = sorted(realc)
        pts = [0.5 * (realc[i] + realc[i + 1]) for i in range(len(realc) - 1)]
        signs = [(-1) ** n * lead_sign]
        for p_ in pts:
            s_ = H.sign_at(p_)
            if s_ == 0:
                undecided = True
                break
            signs.append(s_)
        signs.append(lead_sign)
        nsc = sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])
        if nsc == n and not undecided:
            res = dict(status="all-real", n_real_lb=n, n_nonreal_lb=0, exact=True, prec=prec)
            return res
        boxes = []
        rad = arb(2) ** (-(int(0.8 * n) + 60))
        for c in sorted(cplx, key=lambda v: abs(v.imag)):
            cc = H.newton(complex(c))
            if not (cc.imag.mid() > 0):
                continue
            if any(abs(cc - b).mid() < 1e-6 for b in boxes):
                continue
            if H.krawczyk(cc, rad):
                boxes.append(cc)
                if not want_full:
                    break
        if boxes:
            ok = all(abs(boxes[i] - boxes[j]).mid() > 4 * float(rad.mid()) for i in range(len(boxes)) for j in range(i))
            nrl = nsc if not undecided else 0
            res = dict(status="nonreal", n_real_lb=nrl, n_nonreal_lb=2 * len(boxes) if ok else 2,
                       exact=ok and (nrl + 2 * len(boxes) == n), boxes=boxes, prec=prec)
            if res["exact"] or not want_full:
                return res
        prec *= 2
    return res


def _critical_point(H, x0, iters=60):
    x = arb(x0)
    for _ in range(iters):
        Q, D1, D2 = H.eval(arb(x.mid()), deriv=2)
        step = (D1 / D2).mid()
        x = arb(x.mid()) - step
        if abs(step) < arb(2) ** (-(H.prec - 40)) * (1 + abs(x.mid())):
            break
    return float(x.mid())
