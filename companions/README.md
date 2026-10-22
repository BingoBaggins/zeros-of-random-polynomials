# Companion papers to annals/main.tex — assembled 22 October 2026

Three papers holding the material cut from the five drafts in `../paper/` when `../annals/main.tex` ([M]) was formed
(see `../paper/ANNALS-PLAN.md`, Section 5). Author B. Baggins. Each builds with

    cd companions && tectonic -k --keep-logs paper-X.tex

Build all three twice: they cite each other and [M] through xr-hyper (`\externaldocument[M-]{../annals/main}`,
`[A-]{paper-A}`, …), so the aux files must be current. State at assembly: no warnings, no undefined references,
no `??` in the text; paper B has one harmless underfull box on page 22.

| Paper | File | Pages | Contents (source sections) |
|---|---|---|---|
| [A] Explicit onsets, rates and the tail of the collision constant | `paper-A.tex` | 25 | Paper I §6.1 (explicit onset); Paper II §2.6 (rate n^{-α}, overcrowding, distortion); Paper III §3 (explicit constants, two-point functions, modulus gap, tail of D); questions q:distortion, q:Dtail |
| [B] Nonreal zeros near the collapse time | `paper-B.tex` | 23 | Paper IV as it stood, citations redirected to [M], [A], [C]; Hermite zeros renamed x_k to avoid the clash with the coefficients ξ_k |
| [C] Exact counts beyond the Weyl weight | `paper-C.tex` | 35 | Paper I §6.2 (symbols depending on n), §7 and Appendix A (varying weights); Paper V (profiles, elliptic polynomials, Laguerre flow); I's "Where the Hermite structure is used" and q:varying |

`NOTES-A/B/C.md` record, per paper, what was written, which sentences were reworded, notation aligned with the
22 October changes in [M] (κ→ϰ, R→ρ_q), bibliography entries removed, labels referenced across papers, and the
remaining `TODO(author)` items. Labels keep their original names, so the sources can be compared line by line.

Cross-paper labels that must survive edits: [C] and [B] cite [A]'s `thm:quant`, `thm:quant-main`, `thm:rate`,
`def:condition`, `lem:outer`, `lem:bulk-explicit`, `eq:bulk-explicit`.

Not changed: `../annals/main.tex` still lists [A], [B], [C] as "in preparation"; the companions cite each other
and [M] as "preprint (2026)".
