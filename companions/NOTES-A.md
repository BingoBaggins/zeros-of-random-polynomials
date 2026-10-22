# Notes on paper-A.tex (companion [A])

Status (22 Oct 2026): builds cleanly with `tectonic -k --keep-logs paper-A.tex`. **25 pages**, no warnings, no undefined references, no overfull or underfull boxes, and `pdftotext paper-A.pdf - | grep -c '??'` gives 0.

## What was written
- **Abstract** (new). It covers the explicit onset, the rate n^{-alpha} for alpha < 3/50, the explicit constants of the collision theorem, and the tails of D (exponent 1/900 on {T_* >= 1+eps}, and exp(-c(log x/log log x)^{1/2}) without that restriction). It contains no citations.
- **Introduction**, with these subsections:
  - 1.1 Setting: Weyl polynomials, the heat flow, Lambda_n^std, q, R_*, T_*. One paragraph each on [M, Theorem 1.1] (`M-cor:weyl`) and [M, Theorem 1.4] (`M-thm:main`), plus the identity [M, (1.6)] (`M-eq:heat-hermite`) and the reversed coupling P_{n,t}, Lambda_n ([M, Theorem 1.5], `M-thm:heat`).
  - 1.2 Theorem 1.1 (`thm:quant-intro`): the Paper I statement, unchanged, followed by the paragraph "This is Theorem 2.6 ...". I added one clause, "and they need not be simple", which the hypotheses of `thm:quant` support (multiplicities are counted, and nu_d > 0 excludes only multiple real zeros). [M, Remark rem:hypotheses] says the same.
  - 1.3 Theorem 1.2 (`thm:rate-intro`): the Paper II statement, unchanged, followed by the pointer "part of Theorem 3.6".
  - 1.4 A recap in words of [M, Theorem 9.1] and the definition of D, then Theorem 1.3 (`thm:explicit-intro`, a new label). It restates in words (a) `thm:collision-explicit`, (b) `thm:Dtail` and (c) `cor:Dtail`, with the exponents unchanged: 6eps^{-1}Pi, Pi^{5/2}, Pi^{10}, beta_eps >= eps/(35640 log(40/eps)), 1/900, and (log x/log log x)^{1/2}. The pointer sentence follows: "Part (a) is Theorem 4.1, part (b) is Theorem 4.12, and part (c) is Corollary 4.13."
  - 1.5 Ideas of the proofs. The paragraph on Section 4 adapts Paper III's "Ideas of the proofs".
  - 1.6 Related work, kept brief: it refers to the introduction of [M] for the literature and cites [B], [C], [Ka], [MS], [Fo] and [AW].
  - 1.7 Organisation.
- The commented block "MATERIAL FOR THE INTRODUCTION" is deleted.

## Reworded sentences in the body
- Section 2, first sentence: "The proof above obtains $n_0$ by compactness. ... instead of disks around the $\zeta_i$" is now "The proof of [M, Theorem 1.4] in [M, Section 6] obtains n_0 by compactness. ... we control the deformation a_s(w)=p_in(w)p_out(sw) of that proof ... instead of the disks about the points zeta_i=w_i+w_i^{-1} used there".
- Section 2, "Size of the onset": "has a polynomial small-ball bound Lemma~\ref{lem:nu-small}" now has the lemma in parentheses.
- Section 3, first paragraph: "gives a rate in [M, Thm 1.1]" is now "gives a rate of convergence in ...". I added a pointer to the symbol a^{(t)} of [M, (1.6)] and the sentence "As in [M], we work with the polynomials P_{n,t} of [M, (1.5)] ... real-rooting time Lambda_n ([M, Theorem 1.5])". "In this subsection $q$ is the Gaussian power series" is now "In this section q is the Gaussian power series [M, (1.3)]".
- Remark `rem:rate`: "the only place where an earlier form of the argument lost a factor depending on eps was the bound ..., which gave an exponent of order ..." is now "the only place where the argument could lose a factor depending on eps is the factor Upsilon^k, and its bound ... would give an exponent of order ...". This removes the reference to an earlier draft.
- Proposition `prop:distortion`, first line: "be the distortion of Lemma 2.2(c) for (C_0,R_0)=..." is now "be the distortion of the symbol a^{(t)}, defined in Lemma 2.2(c) with (C_0,R_0)=...". This was needed to fix a 4.8pt overfull box.
- Section 4, opening: added "The notation is that of [M, Section 9]: q ..., P_{n,t}, Lambda_n ..., u_*=R_*e^{i theta_*} ...".
- Section 4.1: "q satisfies (G1)--(G3) with $R\ge1$" is now "q satisfies the conditions (G1)--(G3) of [M, Section 9], with radius of convergence rho_q >= 1". "recall gamma_*=u_*q'(u_*)" is now "... from [M, Section 9.3]".
- Proof of Theorem 4.1, Step 1: "this is (C2). (C1) holds because" is now "this is condition (C2) of [M, Section 9.3]. Condition (C1) there holds because".
- **Flagged spot.** The bare reference was in `\cite[Sections~\ref{M-subsec:collision-rate}--\ref{M-subsec:collision-proof}]{Main}`, in the first sentence of the proof of Theorem 4.1. It now reads "We go through the proof of [M, Theorem 9.1], which runs from [M, Section 9.2] to [M, Section 9.5], and bound each constant." This also removed the build's only warning. With xr-hyper and hyperref under xetex, two external `\ref`s inside one `\cite[...]` raise "hyperref Warning: Invalid page number (0) for remote PDF file": `\Hy@MakeRemoteAction` sets `\Hy@href@page` to 0 without a group. **Other agents: avoid two `M-` refs inside one `\cite[...]`.**
- FIXME: none found. Self-citations (`CompA`): none.

## Notation aligned with [M] (checked against annals/main.tex)
- Expected number of zeros of q in the upper half-disk: kappa(r) is now \varkappa(r). This covers 7 occurrences in Section 4.2 (proof of Theorem `thm:gap` and the paragraph and display before `eq:Rsmall`) and Section 4.3 (proof of `lem:collision-smallball`(a)).
- Radius of convergence of q: R is now rho_q. Section 4.1 setting: "(G1)--(G3) with rho_q >= 1" and "zeros of q in {|u|<rho_q}". Step 1: "1 <= rho_q". Proof of `thm:Dtail`: "hold with rho_q = 1".
- The other renamed symbols do not occur in the text of [A]: varpi (now vartheta), the quadrature nodes xi_i (now x_i), kappa_* (now c_{12}), and the moduli G_{mn} and G(m). G_{mn} occurred only in the deleted commented block. G_j in the proof of `lem:bulk-explicit` is [M]'s transfer matrix, which is unchanged in [M] (line 501).
- For the coordinator: in the current main.aux, `cor:collision-weyl` is Corollary **9.16** (not 9.14; 9.14 is `cor:collision-gauss`), and `thm:sharp` is Theorem **7.13**. The references resolve automatically.

## Bibliography
- Kept: Main [M], CompB [B], CompC [C], An, AW, Fo, Fuj, Ka, LPR, MS.
- Removed (not cited): BI, Bi, Bi2, Bl, Bou, CSV, dB, De, DKMVZ, DLMF, DS, Du1, Du2, Fe, GSC, HHJK, HHJK2, HKPV, KMVV, KN, KVA, Le, LMS, Ma, Me, Ne, Ob, OT, Pe, PM15, PV, Ra, RT, RW, SaTo, Si, SNWW, ST, Sz, vN, Wat, We.
- The `% TODO(author)` on Bl (CRM publication data) was deleted with the entry. No kept entry carried a TODO(author) comment in the bibliography.

## References into [B] and [C]
- No labelled references. [B] and [C] are cited only as `\cite{CompB}` (times t_n -> 1, Section 1.6) and `\cite{CompC}` (symbols depending on n, Sections 1.2 and 1.6), so nothing depends on their aux files.
- In the other direction, paper-C.tex cites these labels of [A], all of which are kept: `A-thm:quant`, `A-def:condition`, `A-lem:outer`, `A-lem:bulk-explicit`, `A-eq:bulk-explicit`.

## Remaining TODO(author) items (all kept)
1. Line 43: address and e-mail.
2. Proof of `lem:annulus`: page-check [AW, Theorem 6.2] in dimension one, and the page in [Ka] for the density 1/(pi(1-x^2)).
3. Section 4.2 (Rice formulas): check the pinpoints Theorems 6.2 and 6.4 of [AW]. Only [LPR, App. A.1, Thm 1] was read.
4. Proof of `lem:nearcritical`: the pinpoint in [An] (Anderson's inequality) has not been checked.
5. Statement on the use of AI: it must describe a line-by-line verification before submission. The section text itself is unchanged.
