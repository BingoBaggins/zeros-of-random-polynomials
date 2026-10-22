# The real-rooting time of heat-evolved random polynomials

## Files

| file | content |
|---|---|
| `main.tex` | the paper (single source file, bibliography included) |
| `main.pdf` | the compiled paper, 67 pages |
| `main.aux`, `main.log`, `main.out` | build by-products |

## Building

Run inside this folder:

    tectonic -X compile --keep-logs --keep-intermediates main.tex

Tectonic leaves a `main.xdv` behind; delete it. The last build (22 October 2026) has no warnings, no overfull or
underfull boxes and no undefined references.

## Structure of `main.tex`

| lines | content |
|---|---|
| 57 | §1 Introduction |
| 274 | Part I, Exact counts of nonreal zeros: §§2–6 (Preliminaries; backward recurrence; bulk; exterior; proof of the exact count) |
| 687 | Part II, The limit law: §7 Heat flow and random power series, §8 The standard coupling |
| 1440 | Part III, The final collision: §9 The final collision, §10 The first-order law |
| 2501 | §11 Remarks and questions |
| 2555 | Statement on the use of AI, then the bibliography |

The paper was assembled from earlier drafts. Four `% MERGE:` comments mark the sentences rewritten for the merge.

## Open items

`main.tex` has nine `% TODO(author)` comments:

- the author's address and e-mail (line 43);
- a comparison with the proof in [Pe] (line 260);
- unchecked citation pinpoints in [Sz], [Ob] and [Si] (lines 306, 715, 1042) and in the bibliography (line 2567);
- where the verification suite is deposited (lines 2526 and 2563);
- the Statement on the use of AI, which must describe the author's line-by-line verification of every proof
  (line 2557).
