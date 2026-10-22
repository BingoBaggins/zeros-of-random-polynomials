# Verification suite for `main/main.tex`

Mechanical checks of the identities and explicit constants, and numerical tests of the sharp predictions, of
*The real-rooting time of heat-evolved random polynomials* (`main/main.tex`, 67 pp.). The scripts are meant to be
deposited with the paper (arXiv ancillary files).
The scripts do not modify the paper.

**Result: every check passes; no discrepancy with the paper was found.** Every number the paper prints that these
scripts recompute is reproduced to the printed digits. A few observations that are not discrepancies are listed under
[Observations](#observations). The table of all checks, with their numbers, is generated in `results/SUMMARY.md`.

## Running

```bash
cd verification
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python run_all.py
```

`run_all.py` runs the six scripts below, writes a log and a JSON file per script to `results/`, and writes
`results/SUMMARY.md`. It exits nonzero if any check fails. A full run takes about 9 minutes on a laptop
(`--quick` is a smoke test with smaller sizes, about 2 minutes; in it the asymptotic checks of Parts 3a and 3b are reported but not judged). Tested with Python 3.14, SymPy 1.14, mpmath 1.3,
python-flint 0.9 (Arb ball arithmetic), NumPy 2.5 and SciPy 1.18. All random inputs use fixed seeds.

| script | content |
|---|---|
| `part1_symbolic.py` | Part 1: symbolic identities (SymPy), identities at random 50-digit points (mpmath), and random tests of explicit inequalities |
| `part2_constants.py` | Part 2: the onset `n_s` of `eq:std-onset` and the other constants of `prop:early` (Arb and exact rationals) |
| `part2_fredholm.py` + `fredholm_lib.py` | Part 2: the Fredholm determinant `P(T_* <= t) = det(I - X(r))`: certified Arb enclosures, Bornemann quadrature, double precision up to r = 0.99, Monte Carlo |
| `part3a_exact_count.py` + `certify_lib.py` | Part 3a: `thm:main` for three admissible symbols, with certified counts of real and nonreal zeros |
| `part3b_firstorder.py` | Part 3b: `Lambda_n` by certified bisection in `t`, against `thm:firstorder-intro`, for the benchmark and two other symbols |
| `part4_exact_checks.py` | Part 4: the exact rational-arithmetic computations of `rem:bernoulli-general` (Remark 11.3), ported from `exact_checks.py` |
| `common.py`, `run_all.py` | reporting and the driver |

**What PASS means.** Each check states its method:
- *symbolic*: a SymPy identity, proved for all values of the free symbols, for all integer parameters up to the
  stated size (for example `n <= 9`);
- *exact rational arithmetic*: a statement decided exactly over Q (for example Sturm counts);
- *interval arithmetic (Arb)*: a rigorous enclosure; a PASS is a proof of the stated numerical fact, modulo the
  correctness of Arb and of the short reductions documented in the code;
- *high-precision random points*: an identity checked at random points to 25–50 digits. This is evidence, not a proof;
- *numerical* / *Monte Carlo*: floating-point evidence. An inequality with an explicit constant that is tested at
  random points can be refuted by such a test, but not proved.

## Part 1 — symbolic identities (`part1_symbolic.py`, 63 checks, all PASS)

Last full run: **140 checks passed, 0 failed** (plus 17 informational lines), `results/SUMMARY.md`.

| claim | label in `main.tex` | method | result |
|---|---|---|---|
| `He_k = exp(-d^2/2) x^k` and `exp(-t/(2n) d_z^2)(sqrt n z)^k/sqrt(k!) = t^(k/2) h_k(sqrt(n/t) z)` | `lem:heat`(a) | symbolic in `n,t,z`, `k <= 14` | PASS |
| heat–Hermite identity `P_{n,t} = t^(n/2) sum_j eta_j t^(-j/2) h_{n-j}(sqrt(n/t) z)` | `eq:heat-hermite`, `eq:PQ` | symbolic in `eta_j,t,z`, `n <= 9` | PASS |
| removal of leading zero coefficients | `eq:reduction` | symbolic, three cases | PASS |
| `exp(-a d^2)` keeps degree and leading coefficient and does not increase the number of nonreal zeros | `lem:heat`(b) | exact Sturm counts, 300 random integer polynomials and rational `a` | PASS |
| Christoffel–Darboux formula | `lem:cd` | symbolic, `n <= 12` | PASS |
| continuant formulas for `A_{n,j}`, `B_{n,j}` | `lem:continuant` | symbolic, `j <= 12`, with `sqrt(beta_r)` replaced by free symbols (a stronger identity) | PASS |
| `Pi_m(y;1,...,1) = U_m(y/2)` (limits `A_j`, `B_j`) | `eq:fixed-j` | symbolic, `m <= 12` | PASS |
| `h_{n-j} = A_{n,j} h_n + B_{n,j} h_{n-1}` | `eq:AB` | random complex points, `n` in {7, 18, 30}, all `j <= n` | PASS |
| defect (matching) expansion | `lem:defect` | symbolic in `delta_e`, `y`, `m <= 9` | PASS |
| first-order coefficients `A^(1)_j`, `B^(1)_j` | `lem:secondorder`(b) | series in `1/n`, symbolic, `j <= 9` | PASS |
| recurrence for `Y_j` | `lem:secondorder`(d) | symbolic in `w`, `j <= 12` | PASS |
| bounds (a), (b), (c) of `lem:secondorder`; `lem:chebyshev`; `lem:transfer`; `lem:head-rate`(a),(b) | same | random points (max ratio lhs/rhs reported, all < 1) | PASS |
| particular solution `P(j)`, `alpha_1`, `alpha_2`, and `F_1^* = gamma(1/(4 omega^2) - (i/4) cot(theta) rho)` | `lem:envelope`(b) | symbolic; plus the series `sum_j a_j Y_j` for a concrete `q` to 30 digits | PASS |
| `F_0(T*,y*) = 0`, `d_t F_0 = -gamma/(2T*)`, `d_y F_0 = beta` | `lem:envelope`(b) | random/high precision, concrete `q` | PASS |
| window identity `beta(zeta - i lambda)/omega = beta xi - gamma tau/(2T*) + F_1^*` | `lem:envelope`(c), `eq:window` | symbolic | PASS |
| `A_j + w B_j = w^j`, `H_a = [[A,B],[-B,A+zB]]`, `det H_a(z) = a(w) a(1/w)` | `lem:limit` | symbolic for polynomial `a` of degree `<= 8`; random points for a power series | PASS |
| colleague matrix `det(zI - C) = sqrt(n!)/(a_0 n^(n/2)) Q(sqrt n z)` | `lem:colleague` | symbolic `n <= 5`; random points `n = 12, 25` | PASS |
| model identity `Phi_c = R(s)[(zeta,lambda) - gamma(s)]`, `|(f, d f)|^2 = (J(s)-c)^2 + (lambda + cos^2 s)^2` | `lem:model`(a) | symbolic | PASS |
| values at the threshold: `f'' = -2 sin s_c`, `f''' = -2 cos s_c`, `det D Phi = 2 sin^2 s_c` | `lem:model`(b) | symbolic | PASS |
| Kepler threshold: `Psi - sin Psi = 2 dist(c, pi Z)`, `cos 2s_c = cos Psi`, `sin 2s_c = sgn(c~) sin Psi`; Dini endpoints `lambda_0 = -1`, `lambda_{pi/2} = 0` | `eq:kepler`, `eq:kepler-intro`, `lem:model` | random points; symbolic for `J`, the cycloid and `1 + 2 lambda_c = -cos 2 s_c` | PASS |
| all zeros of `f_{lambda,c}` real iff `lambda >= lambda_c` | `lem:model`(c) | argument principle and sign changes, 12 random `c`, `lambda = lambda_c +- 0.05` | PASS |
| constants of `lem:model`(d) (`C_J < 2.22`, 11.7, 7.1) and the localisation bounds | `lem:model`(d) | Arb for the constants; 4000 random points for the bounds | PASS |
| model is a backward heat trajectory | `lem:model` (remark after it) | symbolic | PASS |
| `|Phi_c| >= 1/2` on the boundary of `B_0`; arithmetic of the double-zero bound | `prop:threshold`, `thm:firstorder`(c) | random points; exact | PASS |
| heat polynomials `H^+_k = He_k`, `H^-_k(v) = i^-k He_k(iv)`, local scaling at a multiple zero; Jacobian `-P''^2/(2n)` | `lem:localform`, `rem:localform` | symbolic, `k <= 14` | PASS |
| `varkappa(r)`, its derivative and integrated form; `pi varkappa = 4/3 (r^3+r^5) + 152/105 (r^7+r^9) + O(r^11)`; intensity = Matsumoto–Shirai at `n = 1`; angular integral `(2/rho) arctan rho` | `lem:pairs`, `prop:law` | symbolic (series, derivative); 40-digit quadrature | PASS |
| tail constants `4/(3 pi)` (Gaussian), `2/9` (uniform); area `8c^2/3` of `D_0(c)`; box (eq:box) | `thm:tail` | symbolic | PASS |
| joint-law density `(pi W)^-1 ((W+s)/(W-s))^(1/2)`, mass 1, mean `W/2`; density `(1 - cos Psi)/pi` of `Psi(U)` | `thm:jointlaw`, `thm:jointlaw-intro` | symbolic; Monte Carlo (20 000 samples, KS distance 0.0067); equidistribution of `tau_c - W cos Psi(c_n^nat)` over `n <= 2e5` for the benchmark (KS distance 3e-5, mean `tau_c + W/2` to 5 digits) | PASS |
| `tau_c`, `W` in both forms; heat-variable conversions; exact formula for `n(z_n - z*)`; interval length `1/(Im u*)^2` | Notation of §10, `thm:firstorder`(b), `thm:firstorder-intro` | symbolic | PASS |
| `g'(m) = phi_m`; Szegő's phase; `c_n^nat` (introduction) `= c_inf - n pi nu` (`prop:phase`) | `prop:phase` | symbolic | PASS |
| `arg c_n(x*)` and `|c_n(x*)|` asymptotics are `O(1/n)` | `prop:phase`, `eq:PR` | numerical, `n` up to 12800 (`n x error` < 0.55) | PASS |
| carrier identities, Hermite ODE, `eq:carrier-Z`, `epsilon(xi)` | `lem:carrier` | symbolic, `n <= 8` | PASS |
| Mehler's formula at complex arguments | `lem:mehler` | random points | PASS |
| benchmark: `u*`, `q'(u*)`, `rho*`, `tau_c`, `W`, `tau_c + W = a/2`, liminf; `sqrt(n!) Q_n^(t) = (x + b sqrt(n/t)) He_{n-1} - chi He_{n-2} = det(x - J')`; `n(a sqrt(n/(n-1)) - a) = a/2 + O(1/n)`; `nu = 0.2215...`, `tau_c = -0.2590...`, `W = 1.2590...` | `rem:benchmark` | symbolic (`n <= 9` for the determinant); numerical | PASS |
| algebra of `n_ext`, `sum (i+1)^2 x^(i+1)`, `|J(w) - 2| = |w-1|^2/|w|` | `lem:ext-quant` | symbolic | PASS |
| `n! <= e sqrt(n)(n/e)^n`, `(n^n/n!)^(1/2) >= e^(n/2)/(e^2 n)^(1/4)` | `rem:full` | high precision, `n < 400` | PASS |

## Part 2 — rigorous constants (`part2_constants.py`, 14 checks; `part2_fredholm.py`, 29 checks; all PASS)

### The onset `n_s` and the constants of `prop:early`

`n_s'` is computed exactly (the condition is `(s/(1-s))^2 s^n <= 1/16` over Q), then `f_s(n_s - 1) < 0 < f_s(n_s)`
is certified with Arb; together with `f_s(1) < 0` and the convexity `f_s'' = (1 + 2 c_12)/(2(n+1)^2) > 0` (checked
symbolically) this identifies `n_s = min{n >= n_s' : f_s(n) > 0}`.

| `s` | `n_s'` (paper) | `n_s'` (computed) | `n_s` (paper) | `n_s` (computed, certified) |
|---|---|---|---|---|
| 3/10 | 1 | 1 | 13 707 | 13 707 |
| 1/2 | 4 | 4 | 25 051 | 25 051 |
| 4/5 | 25 | 25 | 93 949 | 93 949 |
| 19/20 | 169 | 169 | 493 795 | 493 795 |

Also certified or checked symbolically: `f_s(n)/n -> gamma_s(3/4 - c_12(72(1+s)+1/4)) > gamma_s/2`; the identity
for `beta_s^2/(1-s^2)` and `gamma_s <= 1/(256 e)`; `(4e^-3)^(1/2) <= e^-0.8` (0.44626 against 0.44933); the
chi-squared Chernoff bound; the exponent identities of Steps 1–2 (`eq:std-Vinf`); the 12-point average `eq:std-gain`;
`Phi_max <= 36(1+s) gamma_s`; `c_12 = 2/4095`; and that the right side of `eq:std-early` is below 1 only for
`n > 4 log(13)/gamma_s = 27 444.2` at `s = 1/2` ("about 2.7e4"). For `s = 1/2`: `beta_s^2 = 0.0358889`,
`gamma_s = 3.73842e-4`, `c_s = 9.34606e-5`, `C_s = 145.52`.

### The Fredholm determinant (`prop:fredholm`, `rem:fredholm-num`)

Three independent methods (see the docstring of `fredholm_lib.py`):

1. **Certified (Arb).** The entries `G_mn(r)` are computed from an exact series, `G_mn = (2/pi) sum_s T_k(s)
   r^(p+2s)/(p+2s)` with `T_k(s) = sum_{a+b=s} c_a c_b 2k/(k^2 - 4(a-b)^2)`, `c_a = binom(2a,a)/4^a`, from the
   expansion of `|1 - x e^{2i theta}|^-1` (the series itself is checked against quadrature of `S_k` and against the
   three-term relation of `rem:fredholm-num`). Since `|T_k(s)| <= 2|k|`, the series tail is bounded explicitly. The
   matrix is truncated to even indices `m, m' <= M`; with the bound `|X_mm'| <= (2m+1) r^(m+m'+1)/(1-r^2)` from
   the proof of `prop:fredholm`, the discarded part has trace norm at most `tau_M`, and
   `|det(I-X)/det(I-X_M) - 1| <= exp(||(I-X_M)^-1|| tau_M) - 1`, with the norm of the inverse bounded by its
   Frobenius norm in Arb. This gives:

   | `r` | `t = r^-2` | certified enclosure of `P(T_* <= t)` | `tr X(r) = varkappa(r)` |
   |---|---|---|---|
   | 1/2 | 4 | 0.929249993919685 ± 6e-16 | yes (to 1e-22) |
   | 7/10 | 2.0408 | 0.7258052640705 ± 4e-14 | yes |
   | 4/5 | 1.5625 | 0.482244737187 ± 7e-13 | yes |
   | 9/10 | 1.2346 | [0.10868062123685, 0.10868062123699] | yes |

   So the paper's `P(T_* <= 100/81) = 0.10868062` is correct to all printed digits (the certified value is
   0.108680621237 ± 7e-14).
2. **Bornemann quadrature** (Nyström method with a Gauss–Legendre product rule in polar coordinates) for the operator
   `L` of `eq:Lkernel` restricted to even functions, that is, the kernel `2L(z,w)` on the half disk (same nonzero
   eigenvalues as `X(r)`). At `r = 0.9` the error against the certified value is 2.4e-6, 2.1e-8, 2.0e-10, 1.9e-12 and
   1.4e-13 with 16×32 … 48×96 nodes (exponential convergence, as expected for an analytic kernel); at `r = 0.5` it is
   1e-14 with 12×24 nodes.
3. **Double precision**, the matrix `X_M` with `S_k(rho)` from the Fourier coefficients of `|1 - x e^{i phi}|^-1`
   (FFT) and composite Gauss–Legendre quadrature in `rho` graded towards `r`; it agrees with the certified value at
   `r = 0.9` to 1e-12. At `r = 0.99`: `-log P(T_* <= 10000/9801) = 36.180954832` for `M = 1400` and `M = 1750`
   (paper: 36.180955), with componentwise condition number 75 of the determinant. Also at `r = 0.99`:
   `log det(I - S) = -36.2695` and `sum log(1 - r^(4k+1)) = -42.0901` (paper, `rem:lower`(b): -36.269, -42.090); all
   computed eigenvalues of `X(0.99)` are real, those above 1e-12 lie in (0, 0.9332) (`q:bernoulli`).

Further checks in `part2_fredholm.py`: the bound of `thm:lower` at `t = 100/81` is 0.0112148 (paper: 0.0112) and
`P(T_* <= t) >= prod(1 - t^(-2k-1/2))` at the four certified radii; `lem:qform` (the quadratic form identity, to
2e-13, and `0 <= sigma_r <= 1`); Step 1 of `thm:lower` (`S_n <= R_n`); `lem:even-bergman` (`||L - K||_2` = 0.21,
0.59, 0.85 at `r` = 0.5, 0.8, 0.9, against the bound `1 + log(1/(1-r))`), its strip substitution and the integral
`2 - pi/2`; `int_0^oo -log(1 - e^-u) du = pi^2/6`; `lem:log-ineq` (154 000 random points, `p = 2..8`);
`rem:explicit` (`(t-1)(varkappa - ||L||_2^2/2) = 0.27961` at `r = 0.8`; paper 0.2796); `eq:Lkernel` and
`eq:traces`; and a Monte Carlo test of `prop:fredholm` (`P(k_r = j)` for `j < 5` and `E k_r`, at `r` = 0.5, 0.7,
0.9, 6000–20000 samples, all within 4 standard errors).

**What would make the non-certified values rigorous.** For `r <= 0.9` the values above are already certified. For
`r = 0.99` (and `t - 1` down to 0.0026, as claimed in `rem:fredholm-num`) one needs:
(i) a bound for the truncation, which the inequality above gives (`tau_M` is about `r^M` times a polynomial, so `M`
about 3000 at `r = 0.99`);
(ii) bounds for the entries: either the series above in ball arithmetic (about 1e10 operations at `r = 0.99`, feasible
in C with Arb but not in this Python suite), or a Gauss–Legendre quadrature error bound for
`int_0^r rho^(m+n+1) S_k(rho) d rho`, which needs a bound for `S_k` on a Bernstein ellipse around each panel
(`S_k(rho)` is analytic for `|rho| < 1`, and `|S_k(rho)| <= pi/(1 - |rho|^2)` there);
(iii) control of the rounding in the determinant, either by an Arb LU factorisation or by the componentwise condition
number (75 at `r = 0.99`) times the entrywise error.
For Bornemann's method itself, Bornemann's error analysis (F. Bornemann, *On the numerical evaluation of Fredholm
determinants*, Math. Comp. 79 (2010) 871–915) bounds the quadrature error exponentially in the number of nodes once the
kernel is shown to extend analytically to a polyellipse with an explicit bound there; for the kernel `2L(z,w)` in
polar coordinates the extension is limited by the zeros of `1 - z conj(w)`, `1 + z conj(w)` and `1 - w^2` (treating
`conj(w)` as an independent variable), at distance comparable to `1 - r^2`. These constants have not been worked out
here.

## Part 3 — experiments

### 3a. `thm:main`: exact counts (`part3a_exact_count.py`, 5 checks, all PASS)

For each degree the zeros of `Q_n = sum_j a_j h_{n-j}` are certified with `certify_lib.py`: real zeros by sign changes
at Arb-evaluated test points between approximate zeros (each sign change gives a real zero in a distinct interval), nonreal
zeros by a Krawczyk test on a small box off the real axis (unique zero in the box, and its conjugate), and the count
is certified when the two numbers add up to `n`. Approximate zeros come from the colleague matrix of `lem:colleague`.

| symbol | radius | nonreal zeros of `a` in D | certified (real + nonreal) zeros of `Q_n` | observed onset |
|---|---|---|---|---|
| S1 = `(1-2w+2w^2)(1+5w/3)/(1-10w/13)` (zeros `(1±i)/2`, `-3/5`) | 1.3 | 2 | `n-2 + 2` for every `n` in {10, 20, 40, 80, 160, 320, 400} | `n >= 10` |
| S2 = `(1-24w/17+20w^2/17)(1+32w/13+20w^2/13)/(1-20w/21)` (zeros `0.6±0.7i`, `-0.8±0.1i`) | 1.05 | 4 | n=10: 6+4; 20: 10+10; 40: 24+16; 80: 62+18; 160: 138+22; 320: 312+8; 400: 394+6; 600: 594+6; **800: 796+4; 1200: 1196+4** | `n >= 800` |
| S3 = `q(5w/6)`, `q` a Gaussian sample (seed 4; the symbol of `eq:heat-hermite` at `t = 1.44`) | 1.2 | 4 | n=10: 10+0; 20: 18+2; 40: 38+2; **80 … 800: n-4 + 4** | `n >= 80` |

The locations are as in (a), (b): at the largest degree, `n x (distance from zeta_i to the nearest zero)` is 1.4–1.6 (S1),
1.6–4.2 (S2), 2.1–2.8 (S3), so the distance is `O(1/n)`, and no other zero lies outside `[-2, 2]` (in `z = x/sqrt n`).
For `n <= 40` the counts are also confirmed by Arb's root isolation of the monomial-basis polynomial.

### 3b. `Lambda_n` and the first-order law (`part3b_firstorder.py`, 16 checks, all PASS)

`Lambda_n` is bracketed by bisection in rational `t`, each step certified as above (`all real` gives `Lambda_n <= t`,
a certified nonreal zero gives `Lambda_n > t`, by `lem:heat`(c)); brackets have width about 5e-12 after 31–34 steps
(`n <= 408`). For `n` in {800, 1600} the bracket `[t_N - 1e-12, t_N + 1e-12]` around the Newton solution `t_N` of
`Q = dQ/dx = 0` is certified instead; `n = 3200` is Newton only. In every certified case the lower end has exactly two
nonreal zeros (certified count) and the upper end none, as in `thm:collision`(c), and the Newton double zero lies in
the bracket.

*The benchmark* `q = 1 + bw + 2w^2`, `theta_* = 1.1` (`rem:benchmark`, `tab:benchmark`), with the exact phase `c_n`:

| `n` | `dist(c_n, pi Z)` paper / here | `n(Lambda_n - T_*)` paper / here | `tau_c - W cos Psi_n` paper / here | `n x difference` paper / here |
|---|---|---|---|---|
| 100 | 0.671 / 0.6714 | 0.45390 / 0.4538974 | 0.45007 / 0.4500746 | 0.38 / 0.382 |
| 101 | 0.024 / 0.0243 | -1.24112 / -1.2411197 | -1.24775 / -1.2477479 | 0.67 / 0.669 |
| 207 | 1.515 / 1.5151 | 1.00171 / 1.0017067 | 0.99805 / 0.9980490 | 0.76 / 0.757 |
| 408 | 0.009 / 0.00946 | -1.37582 / -1.3758216 | -1.37228 / -1.3722778 | -1.45 / -1.446 |

Also reproduced: `n(Lambda_n - 2) = 1.0017 < 1.0036 = n(2 sqrt(n/(n-1)) - 2)` at `n = 207`, the bound
`Lambda_n <= a sqrt(n/(n-1))` at every tested `n`, and the cusp values 1.62 (`n = 110`) and 1.45 (`n = 408`) of
`rem:sharp`.

*Comparison with* `n(Lambda_n - T_*) ~ tau_c - W cos Psi(c_n^nat)` (`thm:firstorder-intro`), `n x difference`:

| `n` | benchmark | rational `q = (1+2w)(1+4w/3+20w^2/9)/(1-5w/6)` | Gaussian sample (seed 106) |
|---|---|---|---|
| 100 | 0.15 | -1.43 | -41.9 |
| 200 | -0.88 | 1.87 | -69.6 |
| 300 | 0.47 | -1.62 | 14.8 |
| 400 | 0.42 | 1.73 | -5.2 |
| 800 | 0.02 | 1.50 | -28.0 |
| 1600 | -1.77 | 1.10 | -25.1 |
| 3200 | 0.58 | 0.32 | -20.9 |

(`tau_c, W, nu` = -0.2590, 1.2590, 0.2215 for the benchmark; -1.4176, 1.3889, 0.7749 for the rational `q`
(`T_* = 20/9`); 1.0886, 2.2996, 0.9388 for the Gaussian sample (`T_* = 1.82537`, `theta_* = 2.460`, next nonreal
modulus 0.850 against `R_* = 0.740`).) The scaled error `n |difference| d_n^(1/3)` stays bounded for `n >= 800`
(0.67, 1.45, 29 for the three symbols), the larger values at single degrees occur where `dist(c_n, pi Z)` is small
(cusps), and `n(Lambda_n - T_*)` stays in `[tau_c - W, tau_c + W]`. The place formula of `thm:firstorder-intro` for
`n(z_n - z_*)` behaves in the same way (`n x difference` at most 0.3, 0.2, 2.5 for `n >= 800`). This is consistent with
the theorem; a finite computation cannot confirm an `O(.)` statement.

## Part 4 — exact computations of Remark 11.3 (`part4_exact_checks.py`, 13 checks, all PASS)

The paper cites a program `exact_checks.py` for the two computer-assisted statements of `rem:bernoulli-general`. Its
parts A and B are ported here into the suite's format (its part C concerned a finite model that is not in the paper).
All points are Gaussian rationals, so every decision is made in exact arithmetic over Q or Q(i). The correlation
functions are those of `eq:ms-kernel`: `rho_n = prod_j 1/(pi |1-z_j^2|) * i^(-n) Pf[K^c(z_i,z_j)]`.

| claim | method | result |
|---|---|---|
| at `(z1,z2,z3) = (-3/20+i/20, 9/10+7i/50, -93/100+13i/100)` the pair terms `a_ij` are positive and `T^2/(4 a12 a13 a23)` is a rational number `> 1` (3.9154066450...), so no Hermitian determinantal kernel | exact | PASS |
| the four-disk polynomial `P(tau)` has rational coefficients `1 - 38.0308 tau + 534.8086 tau^2 - 3288.5002 tau^3 + 7440.7487 tau^4` (printed: 38.031, 534.81, 3288.5, 7440.7) | exact, digits compared | PASS |
| `disc(P) < 0`, `gcd(P, P') = 1`, and `P` has exactly two real roots (Sturm), so one pair of simple nonreal roots, at `0.0887903 +- 0.0017573 i` (printed `0.08879 +- 0.00176 i`) | exact; roots by mpmath | PASS |
| cross-checks not in the original program: all 27 Pfaffians are of antisymmetric matrices and satisfy `Pf^2 = det` (determinant by SymPy over Q(i)); `rho_2(z,z) = 0` | exact | PASS |

The numerical statement that the nonreal zeros persist up to `delta = 0.0275` and are absent from `delta = 0.03` on is
not checked here.

## Observations

None of these contradicts a statement of the paper.

1. **Late onsets in `thm:main`.** For S2 (radius 1.05, one `zeta_i = -2.031 - 0.054i` at distance 0.06 from the bulk)
   the count is correct only from `n = 800` on (among the tested degrees); at `n = 600` there is still one extra
   nonreal pair, at `z = 1.940 + 0.002i` near the right edge of the bulk, and at `n = 160` there are nine extra pairs.
   The theorem only asserts some `n_0`, which depends on the symbol through `R_0`, `mu` and `delta_0`; experiments
   with symbols of radius close to 1 need degrees in the thousands.
2. **Pre-asymptotic regime for a Gaussian sample.** For the Gaussian sample of 3b, `n(Lambda_n - T_*)` differs from
   `tau_c - W cos Psi(c_n^nat)` by about -0.33 to -0.50 for `n` from 100 to 250, that is, by about a fifth of `W`; the
   difference becomes `O(1/n)` (constant about 25) only from `n` about 300 on. So degrees "up to a few hundred" are not
   in the asymptotic regime for every `q`; this agrees with the comment in `rem:sharp` that the size of the remainder
   depends on `q`.
3. **Numerical practice (not about the paper).** Naive ball arithmetic for the three-term recurrence of `h_k` loses
   about `0.8 n` bits in the bulk (the radii follow the recurrence with absolute values, which grows like `e^(0.55 n)`);
   the certified computations therefore use about `1.6 n + 200` bits and Krawczyk boxes of radius about `2^-(0.8n+60)`.

## Not covered

The proofs themselves (only identities, constants and predictions are checked); the probabilistic statements other than
`prop:fredholm` and `thm:jointlaw` (for example `thm:std`, `thm:tail` beyond its constant, `lem:resample`); the
statements of `thm:sharp` and `thm:lower` beyond the identities and inequalities listed; `rem:explicit` beyond
`r = 0.8`; the fits reported in `q:law` and `rem:sharp` beyond the two cusp values; the 128-degree runs up to
`n = 25607` of `rem:benchmark`/`rem:sharp` (here `n <= 3200`); the range of `delta` in `rem:bernoulli-general`.
