# Hexagon submissions: what to enter

Field-by-field entries for the four papers on hexagonmath.org/submit. Every value comes from the paper itself (title,
abstract, `\keywords`, `\subjclass`), cut to Hexagon's limits of 10 tags and 5 MSC2020 codes. In the abstracts the
papers' private macros (`\PP`, `\eps`, `\D`, `\R`, ...) are expanded to standard LaTeX.

## Before you start

- **Order.** Submit the main paper first. Each companion cites it as "B. Baggins, *The real-rooting time of
  heat-evolved random polynomials*, preprint (2026)", and the main paper cites the companions the same way; once
  Hexagon assigns identifiers, the bibliographies can be updated in a revision.
- **ZIP files** are in the repository root and are not tracked by git. Each holds a standalone source: comments
  are stripped, and in the companions every cross-reference to another of the four papers (`\ref{M-...}` and so on,
  which needs the other paper's build files) is replaced by its number, so Hexagon can compile each paper on its own.
  Each was compiled from the ZIP alone and its text checked against the PDF in the repository.
- **Formalization:** none. **AI provenance files:** none (optional).
- **Licence:** the GitHub repository is CC0 1.0; choose CC0 if the form offers it, to match.
- **Visible in every PDF:** the date line reads "Draft, October 2026", and the Statement on the use of AI says no
  human expert has read the proofs. Change these in the sources before submitting if you want them different.

## Main paper — `hexagon-submission.zip`

**Files:** `main.tex` (comments stripped) and `anc/` with the verification suite: 10 Python scripts and the JSON and log results of the last run (140 checks passed). 67 pages.

**Title**

```
The real-rooting time of heat-evolved random polynomials
```

**Abstract**

```
Let $G_{n,t}=\exp(-\frac{t}{2n}\partial_z^2)G_n$ be a random Weyl polynomial of degree $n$ under the backward heat flow, with independent coefficients of a law that is not a point mass and has finite variance. Hall, Ho, Jalowy and Kabluchko showed that the empirical zero distribution collapses onto the real line at $t=1$. We show that the first time at which all zeros of $G_{n,t}$ are real, a finite-degree analogue of the de Bruijn--Newman constant, converges in distribution to $T_*=R_*^{-2}$, where $R_*$ is the smallest modulus of a nonreal zero of the random power series $\sum_j\eta_jw^j$ with coefficients of the same law; the convergence is almost sure under a coupling of the degrees that shares the top coefficients, whereas under the standard coupling the set of limit points of the real-rooting time is almost surely $[1,\infty]$. The proof rests on an exact count of the nonreal zeros of analytic combinations $\sum_{j\le n}a_jh_{n-j}$ of orthonormal Hermite polynomials: if the real power series $a$ has radius of convergence greater than $1$, no zeros on the unit circle and simple zeros in the unit disk, then for all large $n$ the combination has exactly as many nonreal zeros as $a$ has in the disk. Its method expresses the polynomial and its derivative in the bulk through a $2\times2$ matrix applied to $(h_n,h_{n-1})$, built from only $O(\log n)$ steps of the backward three-term recurrence, with the remaining terms controlled by a Christoffel--Darboux bound. For Gaussian coefficients the law of $T_*$ is a Fredholm determinant, $\mathbb{P}(T_*>t)\sim\frac4{3\pi}t^{-3/2}$ as $t\to\infty$, and $(t-1)\log\mathbb{P}(T_*\le t)\to-\pi^2/12$ as $t\downarrow1$, half of the constant for complex coefficients. At the scale $1/n$ the last pair of nonreal zeros collides inside the real bulk, at relative position $\cos\theta_*$, where $R_*e^{\pm i\theta_*}$ are the nonreal zeros of minimal modulus, and the real-rooting time is within $O(1/n)$ of $T_*$. In a window of size $1/n$ the polynomial is, up to a positive factor and $O(1/n)$, the model $\zeta\cos(\zeta+c)+\lambda\sin(\zeta+c)$, whose threshold for real-rootedness is given by Kepler's equation; the phase $c$ rotates with $n$, so $n(\Lambda_n-T_*)$ oscillates over an explicit interval and has an explicit limit law jointly with $T_*$.
```

**Primary subject:** math.PR — Probability

**Cross-list subjects:** math.CA — Classical Analysis and ODEs; math.CV — Complex Variables; math.FA — Functional Analysis (optional: the Fredholm determinant)

**Tags** (10). The paper lists 11 keywords; "Kepler's equation" is left out.

```
heat flow, real zeros, Weyl polynomials, random power series, de Bruijn-Newman constant, Hermite polynomials, linear combinations of orthogonal polynomials, collision of zeros, Pfaffian point processes, Fredholm determinants
```

**MSC2020 codes** (5, in this order). The paper lists 8; 35K05, 47B10 and 60G15 are left out.

| code | topic |
|---|---|
| 26C10 | real polynomials: location of zeros |
| 60G55 | point processes |
| 30B20 | random power series |
| 30C15 | zeros of polynomials and analytic functions |
| 33C45 | orthogonal polynomials (Hermite) |

## Companion A — `hexagon-companion-A.zip`

**Files:** `paper-A.tex` (single file, comments stripped). 25 pages.

**Title**

```
Explicit onsets, rates and the tail of the collision constant for heat-evolved Weyl polynomials
```

**Abstract**

```
Let $h_k$ be the orthonormal Hermite polynomials for the standard Gaussian weight, and let $a(w)=\sum_{j\ge0}a_jw^j$ be a real power series with radius of convergence greater than $1$, $a(0)\ne0$, no zeros on the unit circle and only simple zeros in the unit disk. For all large $n$ the polynomial $\sum_{j=0}^{n}a_jh_{n-j}$ has exactly as many nonreal zeros as $a$ has in the disk, but the known proof gives the onset by compactness. We make it explicit: if $a$ has no zeros in $\{1-d\le|w|\le1\}$, with $d$ below an explicit threshold, the count is exact as soon as $n$ exceeds a polynomial in $1/d$ and in the inverses of $|a(0)|$, $\min_{|w|=1}|a|$, $\min_{|w|=1-d}|a|$ and $\min_{[-1+d,1-d]}\max(|a|,|a'|)$; no separation of the nonreal zeros of $a$ is needed. By an exact identity this applies to the random Weyl polynomials $G_n$ with independent standard Gaussian coefficients under the backward heat flow $G_{n,t}=\exp(-\frac t{2n}\partial_z^2)G_n$, whose real-rooting time converges in distribution to $T_*=R_*^{-2}$, where $R_*$ is the smallest modulus of a nonreal zero of the Gaussian power series $q(w)=\sum_j\eta_jw^j$. We show that the convergence holds with the rate $n^{-\alpha}$ on $[1+\varepsilon,\infty)$ for every $\alpha<\frac3{50}$, an exponent that does not depend on $\varepsilon$; the dependence on $\varepsilon$ is removed by a super-exponential tail for the number of zeros of $q$ in a disk. Under a coupling of the degrees that shares the top coefficients, the real-rooting time $\Lambda_n$ satisfies $|\Lambda_n-T_*|\le K/n$ for $n\ge n_0$, with $K$ and $n_0$ depending on $q$. We bound $K$ and $n_0$ by polynomials in condition parameters of $q$, with absolute degrees and, when $T_*\ge1+\varepsilon$, coefficients polynomial in $1/\varepsilon$. For the Gaussian power series this gives, on the event $\{T_*\ge1+\varepsilon\}$, a power tail with exponent $1/900$ for the largest scaled displacement $D$ of the final collision of zeros, and without restriction the tail bound $\exp(-c(\log x/\log\log x)^{1/2})$; the inputs include closed forms of the two-point correlation functions of the zeros of $q$.
```

**Primary subject:** math.PR — Probability

**Cross-list subjects:** math.CA — Classical Analysis and ODEs; math.CV — Complex Variables

**Tags** (10). These are the paper's 10 keywords.

```
Hermite polynomials, linear combinations of orthogonal polynomials, nonreal zeros, explicit onset, heat flow, Weyl polynomials, random power series, rate of convergence, collision of zeros, correlation functions of zeros
```

**MSC2020 codes** (5, in this order). The paper lists 6; 35K05 (heat equation) is left out.

| code | topic |
|---|---|
| 26C10 | real polynomials: location of zeros |
| 60G55 | point processes |
| 30C15 | zeros of polynomials and analytic functions |
| 33C45 | orthogonal polynomials (Hermite) |
| 60G15 | Gaussian processes |

## Companion B — `hexagon-companion-B.zip`

**Files:** `paper-B.tex` (single file, comments stripped). 23 pages.

**Title**

```
Nonreal zeros of heat-evolved Weyl polynomials near the collapse time
```

**Abstract**

```
Let $G_{n,t}=\exp(-\frac{t}{2n}\partial_z^2)G_n$ be the real Weyl polynomial of degree $n$ with independent standard Gaussian coefficients, evolved under the backward heat flow. Hall, Ho, Jalowy and Kabluchko showed that its empirical zero distribution collapses onto the real line at $t=1$. In a coupling of the degrees that shares the top coefficients $\eta_0,\eta_1,\dots$, the main paper to which this one is a companion shows that for each fixed $t>1$, almost surely for all large $n$, the nonreal zeros of $G_{n,t}$ correspond one to one to the nonreal zeros of the Gaussian power series $q(w)=\sum_j\eta_jw^j$ in the disk $\{|w|<t^{-1/2}\}$. We let $t_n=1+s_nn^{-1/3}$ tend to $1$. If $s_n/(\log n)^{5/4}\to\infty$, then with probability tending to $1$ the polynomial at time $t_n$ has at least as many nonreal zeros as $q$ has in $\{|w|<t_n^{-1/2}\}$, and the excess is the number of its nonreal zeros in an explicit set: a strip about the bulk of $[-2\sqrt{t_n},2\sqrt{t_n}]$, of height of order $(\log n)/n$ away from the edges, together with two disks of radius of order $s_n^{6/5}n^{-2/3}$ about the edges. Moreover, with probability tending to $1$ the strip carries no nonreal zero, so that the exact count is reduced to the two edge regions. The proof compares the polynomial, in Joukowski coordinates, with $q$ on a contour that approaches the unit circle to within a distance of order $(\log n)/n$ in the bulk. It rests on a bound, of independent interest, for the ratio $h_{m-1}(\sqrt m\,\zeta)/h_m(\sqrt m\,\zeta)$ of consecutive orthonormal Hermite polynomials for $\zeta$ near the cut $[-2,2]$, under a threshold on the distance to the cut that is logarithmic and sharp up to a constant factor; it is proved with Olver's error bounds for the Liouville--Green approximation of Weber's equation. For the strip, each Hermite polynomial is split into its two parts recessive at $\mp i\infty$, and a reflection argument replaces the comparison.
```

**Primary subject:** math.PR — Probability

**Cross-list subjects:** math.CA — Classical Analysis and ODEs; math.CV — Complex Variables

**Tags** (8). These are the paper's 8 keywords.

```
heat flow, nonreal zeros, Weyl polynomials, random power series, Hermite polynomials, ratio asymptotics, Liouville-Green approximation, Weber's equation
```

**MSC2020 codes** (5, in this order). The paper lists 7; 26C10 and 35K05 are left out, keeping the three primary codes and the two codes for the Liouville-Green method that the paper rests on.

| code | topic |
|---|---|
| 30C15 | zeros of polynomials and analytic functions |
| 33C45 | orthogonal polynomials (Hermite) |
| 60G15 | Gaussian processes |
| 34E20 | WKB / Liouville-Green approximations |
| 41A60 | asymptotic approximations |

## Companion C — `hexagon-companion-C.zip`

**Files:** `paper-C.tex` (single file, comments stripped). 35 pages.

**Title**

```
Exact counts of nonreal zeros beyond the Weyl weight: varying weights, coefficient profiles and the Laguerre flow
```

**Abstract**

```
Let $h_k$ be the orthonormal Hermite polynomials for the standard Gaussian weight, and let $a(w)=\sum_{j\ge0}a_jw^j$ be a real power series with radius of convergence greater than $1$ and $a(0)\neq0$ that has no zeros on the unit circle and only simple zeros in the unit disk. A companion paper shows that for all large $n$ the polynomial $\sum_{j=0}^{n}a_jh_{n-j}$ has exactly as many nonreal zeros as $a$ has in the open unit disk. We extend this exact count to symbols that depend on $n$, under a tail condition matched to the Christoffel--Darboux bound, and to the orthonormal polynomials of varying weights $e^{-nV}$ with $V$ a one-cut regular polynomial.

Under the backward heat flow $\exp(-\frac t{2n}\partial_z^2)$ the zeros of a random polynomial $\sum_{k\le n}\xi_kc_{k,n}z^k$, with independent coefficients of a law that is not a point mass and has finite variance, become real; the real-rooting time is the first time at which all of them are real. For the Weyl profile $c_{k,n}=n^{k/2}/\sqrt{k!}$ the companion paper shows that it converges in distribution to $T_*=R_*^{-2}$, where $R_*$ is the smallest modulus of a nonreal zero of the random power series $\sum_j\xi_jw^j$. Under the heat flow every real polynomial becomes a combination of Hermite polynomials whose coefficients depend on $n$ through the profile. For other coefficient profiles $\exp(ng(k/n)+h_n(k/n))$, with $g$ concave and $h_n$ of logarithmic size, the collapse time $t_{\rm Wig}$ of Hall, Ho, Jalowy and Kabluchko is exactly the time from which these coefficients decay geometrically, and the real-rooting time converges in distribution to $\max(t_{\rm Wig},\rho^2T_*)$, where $\rho$ is the outer radius of the initial zero distribution; subexponential irregularities of the profile can change this limit. Elliptic polynomials reduce exactly to the Weyl case for a Gaussian entire function. Among the flows $\exp(-\tau A)$ with $Ax^k=c_kx^{k-1}$, the Laguerre flows $\exp(-\tau(x\partial_x^2+(\alpha+1)\partial_x))$ are the only ones that produce orthogonal polynomials. We extend the exact count of nonreal zeros of analytic combinations of Hermite polynomials to Laguerre polynomials, with a hard edge for fixed $\alpha$ and two soft edges for $\alpha$ proportional to $n$. The Laguerre flow preserves the polynomials with all zeros in $[0,\infty)$ for every $\alpha>-1$, and under a coupling of the degrees that shares the top coefficients the corresponding real-rooting time has an explicit almost sure limit, which also sees the negative zeros of the random power series.
```

**Primary subject:** math.CA — Classical Analysis and ODEs

**Cross-list subjects:** math.PR — Probability; math.CV — Complex Variables

**Tags** (10). These are the paper's 10 keywords.

```
Hermite polynomials, Laguerre polynomials, varying weights, Riemann-Hilbert problems, linear combinations of orthogonal polynomials, heat flow, coefficient profiles, real-rooting time, Laguerre flow, hard edge
```

**MSC2020 codes** (5, in this order). The paper lists 7; 35K05 and 60G15 are left out.

| code | topic |
|---|---|
| 26C10 | real polynomials: location of zeros |
| 33C45 | orthogonal polynomials (Hermite, Laguerre) |
| 30C15 | zeros of polynomials and analytic functions |
| 42C05 | orthogonal polynomials, general theory (varying weights) |
| 46L54 | free probability (the Laguerre flow is the finite free heat flow) |
