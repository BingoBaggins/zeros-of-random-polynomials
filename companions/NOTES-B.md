# Notes on companion [B], `paper-B.tex` (22 October 2026)

Build: `tectonic -k --keep-logs paper-B.tex`. Result: 23 pages, no LaTeX/package warnings, no undefined
references, no overfull boxes, `pdftotext paper-B.pdf - | grep -c '??'` = 0. All three aux files
(`../main/main.aux`, `paper-A.aux`, `paper-C.aux`) are present, so there is no missing-aux warning.

One `Underfull \vbox (badness 1571)` remains on page 22: the page ends after Question 5.1 and the
"Statement on the use of AI" heading moves to page 23. It is not caused by these edits (a rebuild with the
original Paper IV introduction gives the same box with the same badness). Removing it would need a manual page
break or a preamble change (`\raggedbottom`), so it was left as it is.

## Sentences reworded (layout of the papers)

1. Abstract. Was: "the companion papers show that for each fixed $t>1$ ...". Now: "the main paper to which
   this one is a companion shows that ...".
2. Section 1.2, before (1.4). Was: "one finds ([M, Lemma 7.1(a)]; see (1.4)) that". The script had turned
   `\cite[\eqref{II-eq:PQ}]{PartII}` into this paper's own `\eqref{eq:PQ}`, so the text pointed at the
   equation that follows it. Now: "([M, Lemma 7.1(a)]; see also [M, (7.1)])". The two references sit in
   separate `\cite`s because two remote links in one optional argument gave hyperref "Invalid page number"
   warnings.
3. Theorem A: the statement is now word for word that of `thm:main` ([M, Theorem 1.4]); only the environment
   differs (`theoremA`). The sentence of [M] that defines $w_1,\dots,w_k$ and $\zeta_i=w_i+w_i^{-1}$ was added
   before it, as [M] has it. The old wording was a paraphrase: "for each zero $w$ of $a$ in $\D$ exactly one of
   them lies within distance $\eps$ of $w+w^{-1}$ ...". Heading kept: "[M, Theorem 1.4]".
4. Paragraph after Theorem A. Was: "For a fixed $t>1$, the companion paper [M, Theorem 1.5] deduces ...". Now:
   "For a fixed $t>1$, [M, Theorem 1.5] deduces ...". Was, at the end: "the final collision of zeros at time
   $\Lambda_n$ is described in [M]". Now: "Under the standard coupling, in which one sequence $(\xi_k)$ serves
   for every degree, the real-rooting time of $G_{n,t}$ converges in law to $R_*^{-2}$ [M, Theorem 1.1]; the
   final collision of zeros at time $\Lambda_n$ and its first-order law are described in [M, Section 9] and
   [M, Section 10]."
5. Related work. Was: "For the literature on linear combinations of orthogonal polynomials we refer to [M],
   and for the real-rooting time, its relation to ... [dB,Ne,RT] and the Pfaffian point process ... [Fo,MS],
   to [M]". The conversion had produced the same citation twice. Now: "For the literature on linear
   combinations of orthogonal polynomials, on the real-rooting time, on its relation to ... and on the Pfaffian
   point process ..., we refer to the main paper [M]". Was: "The companion paper [C] treats real-rooting times
   for other coefficient profiles and for the Laguerre flow." Now: "Of the two other companions of [M], [A]
   makes the onset degree of Theorem A explicit [A, Theorem 2.7], gives a rate of convergence in
   [M, Theorem 1.1] for $t\ge1+\eps$, with constants that are not uniform as $\eps\downarrow0$
   [A, Theorem 3.6], and bounds the tail of the collision constant; [C] treats exact counts for varying weights,
   real-rooting times for other coefficient profiles, and the Laguerre flow."
6. Organisation. Was: "From [M] we use the exterior comparison ... and elementary facts about Hermite
   polynomials; from [M], the identity (1.4) and the first intensity ...". The conversion had produced "[M]"
   twice. Now it is one list: "From [M] we use the exterior comparison behind Theorem A [M, Proposition 5.1],
   elementary facts about Hermite polynomials, the identity (1.4) and the first intensity of the nonreal zeros
   of $q$."

All other converted citations read correctly in their sentences and were left as they are: [M, Lemma 7.1(b)],
[M, Definition 1.3], [M, Proposition 5.1], [M, Section 2], [M, Lemma 2.1] (`lem:zeros`), [M, (2.1)]
(`eq:recurrence`), [M, Lemma 2.4] (`lem:resolvent`), [M, Lemma 2.6] (`lem:constancy`), [M, (5.2)]
(`eq:joukowski-dist`), [M, Lemma 7.9] (`lem:pairs`), [M, Lemma 2.2] (`lem:cd`) and [M, Theorem 1.5] in
Sections 3 and 5. The source contained no `FIXME`. No mathematical statement or proof was
changed.

## Notation

[B] names these quantities of [M] by symbol: $R_*$, $\Lambda_n$, $Q_n^{(t)}$, $a^{(t)}$, $P_{n,t}$, $q$,
$\mathcal J$, $r$, $h_k$. All of them still have these names in [M]. [B] does not use [M]'s
$\varkappa(r)$, $\rho_q$, $\vartheta$ (Section 10), quadrature nodes $x_i$, $c_{12}$, $\mathcal G_{mn}$ or
$\mathfrak g(m)$. Several local symbols of [B] are the same as symbols of [M] but denote other quantities.
They are defined locally and were kept:
- $\varpi_n=\log n+34$ (Sections 3 and 4) is [B]'s own threshold. It is unrelated to the old $\varpi$ of [M]
  Section 10, which [M] now writes $\vartheta$.
- $\varkappa=(1+\frac1{2m})^{1/2}$ (the Weber scaling, Sections 2 and 4) is not [M]'s expected count
  $\varkappa(r)$. [B] never uses that count.
- $\mathfrak g(\tau)=|3\tau^2+2|/(4|\tau^2-1|^{5/2})$ (Lemma 2.3) is not [M]'s $\mathfrak g(m)$.
- $\vartheta$ (the angle in (BT), and the angles $\vartheta_k$ in the proof of Theorem 4.7) is local.
- $\rho_c(u)$ (the intensity in the proof of Theorem 3.1) and $\rho=|w|$ are not [M]'s $\rho_q$.
- $G_j$ (proof of Proposition 3.7) and $\kappa$ (scaling factor in Lemma 3.6, covariance entry in Lemma 3.4)
  are local.
- The zeros of $h_m(\sqrt m\,\cdot)$ (Section 2, (2.1), and the proof of Proposition 3.7) were called
  $\xi_1,\dots,\xi_m$, which clashed with the coefficients $\xi_k$ of (1.1). As the coordinator asked, they
  are now $x_1,\dots,x_m$. This matches [M], which renamed its quadrature nodes from $\xi_i$ to $x_i$.
  Before the change, `grep` found no subscripted $x$ anywhere in paper-B.tex. There were five occurrences:
  two in the sentence that defines them, two in (2.1) and one in the proof of Proposition 3.7. Every
  remaining $\xi$ is a coefficient of the Weyl polynomial. $y_k$ was avoided because $y=|\im\zeta|$ in
  [B]. After the change the build is still 23 pages, with 0 `??` and no warnings; only the underfull box
  described above remains.

## Bibliography

- `CompA` ([A]) is now cited: [A, Theorem 2.7] (`A-thm:quant-main`) and [A, Theorem 3.6] (`A-thm:rate`) in
  Related work. Entry kept.
- No entry was removed. Every remaining entry (M, A, C, AT, dB, DLMF, Fo, HHJK, MS, Ne, Ol, RT, SV, Sz) is
  cited. The `% TODO(author)` comments on DLMF and Ol are kept.

## Labels of [A] and [C] referenced

- [A]: `thm:quant-main` (prints Theorem 2.7) and `thm:rate` (prints Theorem 3.6). If the author of [A]
  renames or removes these labels, [B] gets `??`.
- [C]: none by label. [C] is cited only as `\cite{CompC}`.
- [M]: `thm:main`, `cor:weyl`, `thm:heat`, `def:admissible`, `lem:heat`, `eq:PQ`, `prop:exterior`,
  `sec:prelim`, `lem:zeros`, `eq:recurrence`, `lem:resolvent`, `lem:constancy`, `eq:joukowski-dist`,
  `lem:pairs`, `lem:cd`, `sec:collision`, `sec:firstorder`.

## Remaining `TODO(author)` items (all kept)

1. Address and e-mail (preamble).
2. [Sz, Theorem 8.22.9]: pinpoint not checked (Related work).
3. [Ol, Chapter 6, Theorem 11.1]: cited as [SV, Theorem 1] attributes it; the book was not consulted (after
   Theorem 2.2).
4. [AT, Theorem 2.1.1] (Borell–TIS): pinpoint not checked (Section 3.3).
5. [Sz, Theorem 3.2.2] (Christoffel–Darboux): pinpoint not checked (proof of Lemma 4.6).
6. AI statement: it must describe a line-by-line verification by the author before submission. The section
   itself is unchanged.
7. DLMF: release number and access date (bibliography).
8. Olver: Chapter 6, Theorem 11.1 is cited through [SV] (bibliography).

There is no `TODO(author)` on the status of Theorem 1.3 / `thm:near1-bulk`, in this file or in
`paper/paper-IV.tex`, so none could be kept. **Coordinator:** if such a comment is wanted, its wording has to
come from the author.

## Page count

23 pages.
