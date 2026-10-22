# Notes on paper-C.tex (companion [C])

Title: "Exact counts of nonreal zeros beyond the Weyl weight: varying weights, coefficient profiles and the Laguerre flow".
Build: `tectonic -k --keep-logs paper-C.tex` gives 35 pages, no warnings, no overfull or underfull boxes, no undefined references, and no `??` in `pdftotext` output (22 Oct 2026).

## Structure
- Abstract (new; see below).
- §1 Introduction: 1.1 The Weyl case (definition of admissible, Theorem A, heat-flow background); 1.2 Symbols that depend on n (new); 1.3 Varying weights (new, Theorem 1.1 = `thm:varying-intro`); 1.4 Coefficient profiles (Theorem 1.2 = `thm:profiles-intro`, elliptic polynomials); 1.5 The Laguerre flow (Theorem 1.3 = `thm:lag-count-intro`, Theorem 1.4 = `thm:lag-intro`); 1.6 Related work; 1.7 Organisation.
- §2 Symbols that depend on n (Paper I §6.2): Hypothesis N, Theorem `thm:moving`, Corollary `cor:moving-onset`.
- §3 Varying weights (Paper I §7): Theorem `thm:varying`, Proposition `prop:RC`, Lemmas ladder, outliers, perturbed.
- §4 Coefficient profiles (Paper V §2). §5 The Laguerre flow (Paper V §3).
- §6 Remarks and questions (Paper V), then "Where the Hermite structure is used" (Paper I) and Question `q:varying`.
- Appendix A: top recurrence coefficients of varying weights (Paper I).
- Statement on the use of AI (unchanged), bibliography.

## What I wrote
- **Abstract**: one paragraph built from Paper I's abstract (the definition of the symbol, the exact count attributed to "a companion paper", and Paper I's sentence on symbols that depend on n and varying weights), then Paper V's abstract verbatim except "a companion paper shows" changed to "the companion paper shows".
- **Introduction**:
  - Before Theorem A, the paragraph from [M] §1.2 defining w_i, ζ_i and the Joukowski map, plus "we write r for its inverse" (r was used in §2 and §3 before its old definition in §4.1).
  - Theorem A replaced by the statement of `thm:main` in main.tex word for word; checked by `diff`, which shows only the environment name differs.
  - "This paper asks ..." now also names symbols that depend on n and varying weights.
  - New §1.2 and §1.3 from Paper I §§1.3–1.4. The commented block was deleted. The restatement Theorem 1.1 is followed by "This is Theorem 3.1; for V=x²/2 it is Theorem A (Remark 3.3(a))".
  - Related work: the duplicated "we refer to [M] ... to [M]" became one sentence. I added Paper I's paragraph on DKMVZ/Bleher/Bleher–Its/KMVV, which concerns §3 and Appendix A. Paper I's Peherstorfer and Durán paragraph was not carried over, since it belongs to [M].
  - Organisation rewritten for the §2–§6 + Appendix layout. The "we use" paragraph now reads: one main paper [M] and, in §2 only, [A]. It lists what is used from [M]: Theorem A and the lemmas of [M, §§2–6], the CD bound, the bulk, the exterior, the heat lemma, the fixed series and three facts about q. From [A] it uses Theorem 2.6.
- **Theorem A unification**: every `\cite[Theorem~\ref{M-thm:main}]{Main}` in the body is now `Theorem~\ref{thm:A}`, 7 places: §2 (twice in Step 4, once after Theorem 2.2), §3 opening, Remark RC(a), and §6 (twice). The only remaining direct citation is the heading of Theorem A itself.

## Sentences reworded (originals quoted briefly)
- §1.1: "The companion paper [M] applies Theorem A" → "The main paper [M] applies ...".
- §1.1: "[M, Corollary 1.1]" (label `cor:weyl`) → "[M, Theorem 1.1]".
- §1.4: "replaced by its version for symbols that depend on n Theorem 2.2" → "(Theorem 2.2)". The same fix in §4.2, where it became ", Theorem 2.2".
- §1.5: "weights e^{-nV} with V a polynomial of Section 3" → "of Theorem 1.1". Also "not covered by Section 3" → "not covered by Theorem 1.1".
- §2 opening: "depend on n Lemma profile-identity and Theorem profile." → "(Lemma 4.1 and Theorem 4.2 below)".
- §3 opening: "The proof of [M, Theorem 1.4] uses only a few properties of the Hermite polynomials (Section 6)" → "The proof of Theorem A in [M, Sections 2–6] uses only a few properties ..., which we list in Section 6.5". A new label `subsec:hermite` was added for that subsection.
- §3 proof, Bulk: "Let (a_s) satisfy (H1), (H2), (H4)" → "(H1), (H2) and (H4) of [M, Hypothesis 4.1]". This also removes an overfull box.
- §3 proof, Counting: "which is what (i) and the Rouché step (ii) use" → "... (ii) there use", meaning [M, §6].
- §4.1: heading "Results used from the companion papers" → "Results used". The recall of r, Q_n[b], Hypothesis N (N1)–(N3) and Theorem moving became one pointer paragraph to §2 and §1. The recalls of [M] (fixed power series, heat flow, random coefficients) are kept unchanged.
- §5 opening: added that c, d, ρ_{k,n}, σ_{k,n} and r_V "play the roles they have in Section 3".
- §6.5: "through the coefficient profile Lemma ..." → "(Lemma 4.1)". Also "costs Section 5." → "costs (Section 5)." and "orthogonal polynomials Proposition 5.2." → "(Proposition 5.2)."
- Multiple `\ref`s inside one `\cite[...]` note are now braced, e.g. `{\ref{M-sec:prelim}}--{\ref{M-sec:proof}}`. This fixes the 6 hyperref warnings "Invalid page number (0) for remote PDF file", which come from xr-hyper when two remote links share one group. The headings were not the source of those warnings; all math headings already used `\texorpdfstring`.
- The flagged `M-sec:proof`, `M-lem:transfer` and `M-sec:exterior` lead-ins were checked against main.tex and are correct as they stand: δ_0, ε, D_i, U, K are defined in [M, §6]; C_1 = 4e^{c_κ} matches [M, §4.2]; Lemma 3.1 (continuant) to Lemma 3.3 (transfer); Sections 3–5. No `FIXME` was found.

## Notation alignment with [M] (22 Oct changes)
I checked all of these; none needed a change.
- C names the radius of q as ρ_q in §4.1, which matches [M] (`thm:fixedq`, `eq:Rstar-det`).
- C never uses ϰ(r), ϖ/ϑ, quadrature nodes ξ_i, or the moments 𝒢_{mn}/𝔤(m).
- κ_* in Corollary 2.5 is the κ of [A, Lemma 2.3] at starred parameters. It is not [M]'s old κ_* (now c_{12}).
- G_{n,j}, G_j, H_{n,M}[a], H_a, v_n, c_κ and (H1)–(H4) in §2 and §3 are [M] §4's current symbols.
- The x_i in Lemma 5.9 are local zeros, not quadrature nodes.

## Bibliography
- Removed (uncited): DS, Du1, Du2, GSC, KVA, Pe.
- Kept: [M], [A], [B]. [B] is cited in §1.1 (times close to 1) and §6.5.
- All `% TODO(author)` comments on the remaining entries are kept.

## Labels of [A] and [B] referenced
- [A] (from paper-A.aux, all present): `A-def:condition` (Def. 2.1), `A-lem:outer` (Lemma 2.2), `A-lem:bulk-explicit` (Lemma 2.3), `A-eq:bulk-explicit`, `A-thm:quant` (Theorem 2.6), `A-thm:rate` (Theorem 3.6).
- [B]: no labels, only `\cite{CompB}`.

## Remaining TODO(author) items (14, all kept)
- Author address.
- Br/Ha normalisations.
- SaTo pinpoints (I.1.3, III.2.1, IV.1.10).
- Chi/Sz Favard pinpoints.
- Sz Thm 8.1.3.
- DKMVZ not accessed.
- Lenells pinpoints (checked against arXiv; journal numbering not compared).
- AI statement must describe human verification.
- AKK article number.
- Bleher book data.
- Ch68 last page.
- GM journal data.
- HHJK page range.

## For the coordinator
- The abstract calls [M] "a companion paper" because an abstract cannot cite; the introduction calls it "the main paper [M]".
- Page count rose from 33 to 35 because of the abstract and the new introduction subsections.
