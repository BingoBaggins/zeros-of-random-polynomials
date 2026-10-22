# The real-rooting time of heat-evolved random polynomials

## Overview

Take a random Weyl polynomial `G_n(z) = ∑_{k≤n} ξ_k (√n z)^k/√(k!)` and run it under the backward heat flow
`G_{n,t} = exp(−(t/2n) ∂_z²) G_n`. Hall, Ho, Jalowy and Kabluchko showed that its zeros collapse onto the real line at
time `t = 1`, in the sense of empirical measures. That picture cannot see a bounded number of nonreal zeros, so it says
nothing about the first time `Λ_n` at which *all* zeros are real, a finite-degree analogue of the de Bruijn–Newman
constant. The paper identifies the limit of `Λ_n` and describes how the last pair of nonreal zeros becomes real.

- **Main result (Theorem 1.1).** `Λ_n` converges in distribution to `T_* = R_*^{−2}`, where `R_*` is the smallest
  modulus of a nonreal zero of the random power series `q(w) = ∑_j η_j w^j` with coefficients of the same law. This
  holds for every coefficient law with finite variance that is not a point mass. So the last nonreal pair survives the
  macroscopic collapse by a random amount of time.
- **The engine (Theorem 1.4).** A deterministic exact count: if a real power series `a` has radius of convergence
  greater than 1, no zeros on the unit circle and simple zeros in the unit disk, then for all large `n` the combination
  `∑_{j≤n} a_j h_{n−j}` of orthonormal Hermite polynomials has exactly as many nonreal zeros as `a` has in the disk.
- **The final collision (Theorems 1.6–1.8).** For Gaussian coefficients the last pair collides inside the real bulk,
  within `O(1/n)` of an explicit point, and `n(Λ_n − T_*)` oscillates with `n` over an explicit interval governed by a
  Kepler equation, with an explicit joint limit law.
- **Side results.** For Gaussian coefficients the law of `T_*` is a Fredholm determinant, with
  `(t−1) log P(T_* ≤ t) → −π²/12` as `t ↓ 1`, half of the constant for complex coefficients (Theorem 1.9). Under the
  standard coupling of the degrees `Λ_n` has no almost sure limit (Theorem 1.10).

The identities, the explicit constants and several sharp predictions are checked by a verification suite (140 checks,
all passing). No human expert has read the proofs; see below for exactly what is and is not checked.

## Provenance

The mathematics, the proofs, the verification scripts and the text in this repository were produced almost entirely
by AI sessions (systems from Anthropic, OpenAI and Google). A single person, who prefers to stay anonymous and writes
here as B. Baggins, chose the question, set the scope, directed the sessions, and decided what to keep. No human
expert has read the proofs. The proofs went through several AI referee rounds, and a further AI review of the merged
paper found no mathematical errors and reproduced its numerical claims; none of that is an independent human review
and it should not be treated as one.

The repository is released so that anyone interested can take the ideas, check them, correct them, or carry them
further, with no obligation to the author. If you build on it, a link back to this repository is appreciated but not
required.

## Acknowledgement

The ideas in this repository come from the published work of other people. A language model has no mathematics of
its own: what it knows, it learned from the literature, and what it produced here is a recombination and extension of
results that others found. In that sense the people below are the real creators of this work, and publishing it
anonymously, with no name attached that could take credit from theirs, is the author's way of respecting that.

The paper builds most directly on:
- **Brian Hall, Ching-Wei Ho, Jonas Jalowy and Zakhar Kabluchko**, whose heat-flow model, normalisation and
  macroscopic collapse theorem are the starting point of the whole paper.
- **James Shohat**, whose theorem on quasi-orthogonal polynomials gives the endpoint of the deformation in the proof of
  Theorem 1.4, and **Alan Beardon and Kathy Driver**, **Pál Turán**, **Arieh Iserles and Edward Saff**, **Franz
  Peherstorfer** and **Antonio Durán**, for the theory of zeros of linear combinations of orthogonal polynomials;
  Peherstorfer's exact counts are the closest prior result to Theorem 1.4.
- **Sho Matsumoto and Tomoyuki Shirai**, who found the Pfaffian structure of the zeros of the Gaussian power series,
  **Peter Forrester**, who found it through truncated orthogonal matrices, and **Eric Rains**, whose identity turns it
  into the Fredholm determinant of Theorem 1.9.
- **Yuval Peres and Bálint Virág**, whose hole probability for the complex Gaussian power series is the model for
  Theorem 1.9, and **Mikhail Sodin and Boris Tsirelson**, whose hole-probability argument is behind Proposition 1.11.
- **Larry Shepp and Robert Vanderbei** and **Mark Kac**, for the intensities of the complex and real zeros.
- **Nicolaas de Bruijn, Charles Newman, Brad Rodgers and Terence Tao**, and the **Polymath 15** project, for the
  de Bruijn–Newman constant and the dynamics of its zeros: the local form of a collision (Lemma 9.3) is Polymath 15's,
  and the carrier of the window model is the equilibrium of Rodgers and Tao. **George Csordas, Wayne Smith and Richard
  Varga**, for Lehmer pairs, and **Philip Rosenbloom and David Widder**, for the heat polynomials.
- **W. Shi, Gergő Nemes, Xiang-Sheng Wang and Roderick Wong**, whose error bounds for the Plancherel–Rotach
  asymptotics make the carrier phase of Section 10 explicit.
- **Percy Deift, Thomas Kriecherbauer, Kenneth McLaughlin, Stephanos Venakides and Xin Zhou**, whose Riemann–Hilbert
  asymptotics carry the exact count to varying weights in companion paper C.

The classical tools (the Laguerre–Pólya theory of the heat flow, Christoffel–Darboux and Gaussian quadrature, Mehler's
formula, Rouché's and Hurwitz's theorems, Weyl's equidistribution theorem, the Fredholm theory of trace-class
operators) and the books of Szegő, Watson, Obreschkoff, Simon, Billingsley, Azaïs and Wschebor, and Hough, Krishnapur,
Peres and Virág are older debts of the same kind. The verification suite rests on SymPy, mpmath and Fredrik
Johansson's Arb library (now part of FLINT).

None of these people has seen this repository, contributed to it, or is responsible for anything in it. Any errors are
the author's and the AI's.

## What is checked and what is not

- **Checked mechanically (`verification/`, 140 checks).** The algebraic identities of Sections 2–3 and the window
  model of Section 10, symbolically (SymPy); the onset `n_s` and the constants of Proposition 1.11, and the Fredholm
  determinant values of Remark 7.12, by interval arithmetic (Arb); the exact rational computations of Remark 11.3; and
  the zero counts of Theorem 1.4 for three symbols, by certified root counting. See `verification/README.md` for the
  method behind each check.
- **Tested numerically, not proved.** `Λ_n` computed by certified bisection against the first-order law of
  Theorem 1.7, for three symbols and degrees up to 3200. Agreement with an `O(·)` statement at finitely many degrees is
  evidence, not proof.
- **Not checked by anyone but the AI that wrote them.** All the proofs: the bulk estimate and the deformation behind
  Theorem 1.4 (Sections 3–6), the probabilistic arguments of Sections 7–8, the collision analysis of Section 9 and the
  window expansion of Section 10. A referee would concentrate on Sections 3–6, since everything else rests on
  Theorem 1.4.
- **The companion papers** are not covered by the verification suite.

Treat every statement as a conjecture with a detailed proof until a person with the relevant expertise has read it.

## Contents

| Path | Content |
|---|---|
| `main/main.tex`, `main/main.pdf` | the paper (67 pages); `main/README.md` lists its sections and open items |
| `companions/paper-A.tex`, `.pdf` | [A] explicit onsets, rates and the tail of the collision constant (25 pages) |
| `companions/paper-B.tex`, `.pdf` | [B] nonreal zeros near the collapse time (23 pages) |
| `companions/paper-C.tex`, `.pdf` | [C] exact counts beyond the Weyl weight: varying weights, coefficient profiles and the Laguerre flow (35 pages) |
| `companions/NOTES-*.md` | what was written and changed in each companion paper |
| `verification/` | the verification suite, with the results of its last run in `verification/results/` |
| `LICENSE` | CC0 1.0 Universal |

## The results

`q(w) = ∑_j η_j w^j` has independent coefficients of law `μ` (finite variance, not a point mass), `R_*` is the
smallest modulus of a nonreal zero of `q`, and `T_* = R_*^{−2}`. The reversed polynomials `P_{n,t}` share the top
coefficients `η_0, η_1, …` across degrees; for each `n` they have the law of `G_{n,t}`.
- **Theorem 1.1.** `Λ_n → T_*` in distribution, and `P(G_{n,t} is real-rooted) → P(T_* ≤ t)` for `t > 1`.
- **Proposition 1.2.** `1 < T_* < ∞`; for Gaussian coefficients `P(T_* > t) ~ (4/3π) t^{−3/2}`. The law of `T_*`
  depends on `μ` (Remark 7.10).
- **Theorem 1.4.** The exact count of nonreal zeros of `∑_{j≤n} a_j h_{n−j}`, with one zero near each
  `w_i + w_i^{−1}` and all others real in `(−2−ε, 2+ε)` after scaling.
- **Theorem 1.5.** Under the reversed coupling, `Λ_n → T_*` almost surely.
- **Theorem 1.6.** The final collision: for Gaussian coefficients, `|Λ_n − T_*| ≤ K/n`, and the last multiple zero is
  a double real zero within `K/n` of `2R_*^{−1} cos θ_*`.
- **Theorems 1.7 and 1.8.** `n(Λ_n − T_*) = τ_c − W cos Ψ(c_n) + O(n^{−2/3})`, with `Ψ` the root of a Kepler
  equation; its limit points fill `[τ_c − W, τ_c + W]`, and jointly with `T_*` it has an explicit limit law.
- **Theorem 1.9.** `P(T_* ≤ t) = det(I − X(t^{−1/2}))` for an explicit trace-class `X`, with
  `(t−1) log P(T_* ≤ t) → −π²/12`.
- **Theorem 1.10 and Proposition 1.11.** Under the standard coupling the limit points of `Λ_n^std` are almost surely
  `[1, ∞]`, and real-rootedness before time 1 has exponentially small probability, with an explicit onset.
- **Questions 11.1, 11.2, 11.4 and 11.5.** The second-order term of the hole probability; whether the number of zeros
  in a half disk is a sum of independent Bernoulli variables; ergodic frequencies under the standard coupling; the
  optimal error at the cusp phases.

## Build

**Paper.** Run this in `main/`. `--keep-intermediates` keeps `main.aux`, which the companion papers read:

```bash
tectonic -X compile --keep-intermediates main.tex
```

**Companion papers.** Run this in `companions/` after building the paper. It builds each paper twice, because they
cross-reference each other (see `companions/README.md`):

```bash
for f in A B C A B C; do tectonic -X compile --keep-intermediates paper-$f.tex; done
```

**Verification suite.** Run these in `verification/`. A full run takes about 9 minutes:

```bash
python3 -m venv .venv
```

```bash
.venv/bin/pip install -r requirements.txt
```

```bash
.venv/bin/python run_all.py
```

## Still open

- **An expert read.** The proofs need an expert reader, above all Sections 3–6 (Theorem 1.4), ideally someone working
  on orthogonal polynomials, and Section 7 on the Pfaffian structure, ideally a probabilist working on random analytic
  functions.
- **The author's items.** The `% TODO(author)` comments in `main/main.tex` (listed in `main/README.md`): unchecked
  citation pinpoints, a comparison with Peherstorfer's proof, and where the verification suite is deposited.
- **The questions of Section 11**, and those of the companion papers.

## Licence

Everything in this repository (the papers, the scripts and this text) is released under CC0 1.0 Universal; see
`LICENSE`. You may copy, change, build on and publish any of it, for any purpose, without asking and without
attribution. A link back to this repository is appreciated but not required.
