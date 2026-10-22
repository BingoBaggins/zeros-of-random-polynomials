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
| 56 | §1 Introduction |
| 273 | Part I, Exact counts of nonreal zeros: §§2–6 (Preliminaries; backward recurrence; bulk; exterior; proof of the exact count) |
| 686 | Part II, The limit law: §7 Heat flow and random power series, §8 The standard coupling |
| 1439 | Part III, The final collision: §9 The final collision, §10 The first-order law |
| 2500 | §11 Remarks and questions |
| 2555 | Statement on the use of AI, then the bibliography |

The paper was assembled from earlier drafts. Four `% MERGE:` comments mark the sentences rewritten for the merge.

## Open items

`main.tex` has five `% TODO(author)` comments, all about the literature:

- a comparison with the proof in [Pe] (line 259);
- unchecked citation pinpoints in [Sz], [Ob] and [Si] (lines 305, 714, 1041) and in the bibliography (line 2559).

## Submission

`../hexagon-submission.zip` is the source package for hexagonmath.org: `main.tex` at the top level and the
verification suite in `anc/` (Python scripts and the JSON and log results of the last run; Hexagon accepts only
those file types there). Its `main.tex` has all comments stripped: comment-only lines are dropped and trailing
comments cut to a bare `%`, which leaves the PDF text identical.
