# AGENTS.md

Guidance for coding agents (and people) working in this repository.

## What this is

Quadriceps.jl serves positive-weight cubature rules for two weights: the Gaussian weight (GH,
`ghpos`) and the uniform weight on the cube (Le, `lepos`). The rules are data: 148 text files
under `data/`, selected from the rule bank of the author's designed-quadrature project. The
code is small; it looks rules up, builds tensor products when asked, and checks exactness.

## Layout

| path | what it holds |
|---|---|
| `src/Quadriceps.jl` | module, exports (`ghpos`, `lepos`), `public` names |
| `src/index.jl` | `RuleInfo`, the catalog (`INDEX`, read from `data/index.tsv` in `__init__`), rule-file parser, cache, `available` |
| `src/plan.jl` | which rule answers a request: `beststored`, `atom`, `cheapest` (dynamic program over splits of `d`), the no-rule error, `gauss1d`, `tensor`, `materialize` |
| `src/api.jl` | `ghpos`, `lepos`, `nnodes`, `ruleinfo`: the `(d, q)` methods and the `(d; p)` methods, and the `normalize = false` transforms |
| `src/verify.jl` | `exactness_error`: largest relative monomial error against closed-form moments |
| `data/<family>/*.csv`, `data/index.tsv` | the rules and their catalog — **generated, never edited by hand** |
| `build/build_data.jl` | regenerates `data/` and `docs/src/rules.md` from the rule bank (own environment in `build/`) |
| `docs/` | Documenter.jl site (own environment); `docs/src/rules.md` is generated |
| `test/runtests.jl` | checks every stored rule, both conventions, the fallback, the error paths |
| `NOTICE.md` | third-party credit and the MIT notice that travels with the derived rules |

## Commands

```
julia --project=. -e 'using Pkg; Pkg.test()'                      # about 30 s; checks all 148 rules
julia --project=build build/build_data.jl [sync folder]           # rebuild data/ and docs/src/rules.md
julia --project=docs docs/make.jl                                 # build the docs into docs/build/
```

Use a released Julia (1.10 or later; the author runs 1.13), never a nightly. The build script
needs the project's sync folder (default `~/Dropbox/oldDesignedQuadrature-sync`), which exists
only on the author's machines; everything else works from a plain clone.

## Conventions that must hold

* **`q` and `p`.** The positional argument is `q`, the number of nodes of the one-dimensional
  Gauss rule, as in `gausshermite(q)` and `gausslegendre(q)`; the rule returned has degree
  `p = 2q - 1`, and for `d = 1` it is the `q`-node Gauss rule. The degree is passed as the
  keyword `p` to the one-positional-argument methods (`ghpos(d; p = 7)`). Internally everything
  works in `p`. Do not make `q` mean the number of nodes returned.
* **Normalized frame inside.** Files, cache, tensor products and `exactness_error` all use the
  normalized frame: `N(0, I_d)` for GH, the uniform density on `[0,1]^d` for Le, weights
  summing to 1. `normalize = false` is applied once, at the end, in `src/api.jl`
  (GH: `x/√2`, `w·π^(d/2)`; Le: `2x - 1`, `w·2^d`). One-dimensional GH factors come from
  `gausshermite(q; normalize = true)` — FastGaussQuadrature's default is the other convention.
* **`normalize = true` is the default**, unlike FastGaussQuadrature. This is deliberate.
* **`pragmatic`.** `false`: a request that no stored rule covers throws an `ArgumentError`.
  `true`: the cheapest tensor product of stored rules and one-dimensional Gauss rules. A stored
  rule wins ties. A request is covered by any stored rule of the same `d` and degree `≥ p`.
* **Positive weights only.** Every stored rule has strictly positive weights and a relative
  monomial error below `1e-11`; the build script enforces both and the tests repeat them.
  Never loosen the gate to admit a rule.
* **Returned arrays are copies.** The cache is never handed out.
* **Matrix shape.** `X` is `n × d`, one node per row, also for `d = 1`.

## Data rules

* Do not edit `data/` or `docs/src/rules.md` by hand; change `build/build_data.jl` and rerun it.
  The bank changes over time, so a rebuild may change node counts; update the numbers quoted
  in `README.md` and `docs/src/` when it does.
* A rule is stored only if it has fewer nodes than every tensor product of lower-dimensional
  stored rules and Gauss rules, or attains Möller's bound. The bank's own tensor fills stay out.
* Credit: the `origin` of a rule comes from the project's `rules/literature.tsv` (transcribed,
  same-rule) and `rules/lineage/*.tsv` (derived from Diallo and Worku). A rule that is, or
  descends from, a published rule must never be presented as the author's own. If the set of
  such rules changes, update `NOTICE.md` and `docs/src/credits.md`.

## Style

* American spelling everywhere, identifiers included.
* Unicode operators where Julia accepts them: `≠`, `≤`, `≥`, `∈`, `∉`, `≡`, `≢`.
* Keep the dependency list at FastGaussQuadrature alone; build and docs tools belong in the
  environments under `build/` and `docs/`.
* Docstrings for everything exported or declared `public`; the docs build runs with
  `checkdocs = :exports`.
* When behavior changes, update together: the docstring, `README.md`, `docs/src/guide.md`, and
  the tests.

## Repository

Private, `github.com/NittanyLion/Quadriceps.jl`, branch `main`. No license has been chosen yet;
do not add one without the author. `Manifest.toml` files and `docs/build/` are not committed.
