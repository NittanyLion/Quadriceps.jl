# AGENTS.md

Guidance for coding agents (and people) working in this repository.

## Paper and deposit

The rules are described in Joris Pinkse, *Positive weight Hermite and Legendre quadrature rules* (2026),
on arXiv as **[arXiv:2609.26840](https://arxiv.org/abs/2609.26840)** and deposited on Zenodo. All three Zenodo records are published:

* paper: **[10.5281/zenodo.22904159](https://doi.org/10.5281/zenodo.22904159)** (concept DOI: it always resolves to the latest version)
* rules, to 80 digits: **[10.5281/zenodo.22881864](https://doi.org/10.5281/zenodo.22881864)** (one record for both weight families)
* software snapshot of all five packages, v0.1.0: **[10.5281/zenodo.22883240](https://doi.org/10.5281/zenodo.22883240)**

Cite the paper by its arXiv identifier and its Zenodo DOI together. Keep the block at the top of
`README.md` in step with this one, in all four packages (Quadriceps.jl, quadriceps-py, quadriceps-r, quadriceps-stata).

## What this is

Quadriceps.jl serves positive-weight cubature rules for two weights: the Gaussian weight (GH,
`ghpos`) and the uniform weight on the cube (Le, `lepos`). The rules are data: one binary file,
`data/rules.bin`, selected from the rule bank of the author's designed-quadrature project. The
code is small; it looks rules up, builds tensor products when asked, and checks exactness.

## Layout

| path | what it holds |
|---|---|
| `src/Quadriceps.jl` | module, exports (`ghpos`, `lepos`), `public` names |
| `src/index.jl` | `RuleInfo`, the catalog (`INDEX`, read from `data/index.tsv` in `__init__`), reader and writer of `data/rules.bin` (format `QUADRICEPS1`), cache, `available` |
| `src/plan.jl` | which rule answers a request: `beststored`, `atom`, `cheapest` (dynamic program over splits of `d`), the no-rule error, `gauss1d`, `tensor`, `materialize` |
| `src/api.jl` | `ghpos`, `lepos`, `nnodes`, `ruleinfo`: the `(d, q)` methods and the `(d; p)` methods, and the `normalize = false` transforms |
| `src/verify.jl` | `exactness_error`: largest relative monomial error against closed-form moments |
| `data/rules.bin`, `data/index.tsv` | all rules in one binary file, and their catalog — **generated, never edited by hand**; format in `docs/src/format.md` |
| `build/build_data.jl` | regenerates `data/`, `docs/src/rules.md` and the generated blocks of `README.md` from the rule bank (own environment in `build/`) |
| `build/update.sh` | the unattended updater (hourly cron on one machine): rebuild, test, commit, push, then sync the sister packages `../quadriceps-py`, `../quadriceps-r` and `../quadriceps-stata` |
| `docs/` | Documenter.jl site (own environment); `docs/src/rules.md` is generated |
| `test/runtests.jl` | checks every stored rule, both conventions, the fallback, the error paths |
| `NOTICE.md` | third-party credit and the MIT notice that travels with the derived rules |

## Commands

```
julia --project=. -e 'using Pkg; Pkg.test()'                      # about 30 s; checks all 142 rules
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
  normalized frame: `N(0, I_d)` for GH, the uniform density on `[0,1]ᵈ` for Le, weights
  summing to 1. `normalize = false` is applied once, at the end, in `src/api.jl`
  (GH: `x/√2`, `w·πᵈᐟ²`; Le: `2x - 1`, `w·2ᵈ`). One-dimensional GH factors come from
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

* **One binary file.** The rules are stored in `data/rules.bin` only. Never add per-rule CSV or
  other text files to this repository or to the twins; the history was purged of them on
  2026-09-18 at the author's request. A format change needs a new magic/`fmt`, a matching
  change in all three readers, and `docs/src/format.md`.

* Do not edit `data/`, `docs/src/rules.md` or the `GENERATED` blocks of `README.md` by hand;
  change `build/build_data.jl` and rerun it. The bank changes over time and `build/update.sh`
  commits the changes by itself, so never write a node count of a non-minimal rule into prose:
  put it in a generated block (README) or compute it in an `@eval`/`@repl` block (docs).
* This package is the master copy for the Python, R and Stata twins (`../quadriceps-py`,
  `../quadriceps-r`, `../quadriceps-stata`); their data, `RULES.md` and `NOTICE.md` are copied
  from here. Change behavior in all four together. The Stata twin cannot be tested here (no
  Stata); its data files carry the package's name (`quadriceps_rules.bin`, `quadriceps_index.tsv`)
  because `net install` flattens directories.
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

`authored_by.svg` is the author's shield (the same file as in MemoryLayouts.jl); the README shows
it after the other badges. Do not replace it with a generated shields.io badge.

Quality gate: `Aqua.test_all(Quadriceps)` runs as part of the test suite (ambiguities, stale or
unbounded dependencies, undefined exports, piracy); new dependencies need a `[compat]` entry.
The logo is `docs/src/assets/logo.svg` (Documenter picks it up; the twins carry copies).

GitHub Actions (`.github/workflows/CI.yml`) runs the tests on every push to `main` and on pull
requests; the unattended data updates trigger it too. Check `gh run list` after pushing. The
`docs` job also deploys the Documenter site to the `gh-pages` branch (`deploydocs` in
`docs/make.jl`, pushed with the workflow's own `GITHUB_TOKEN`), which GitHub Pages serves at
`https://NittanyLion.github.io/Quadriceps.jl/dev/`. The branch is generated; never edit it.

Public, `github.com/NittanyLion/Quadriceps.jl`, branch `main`. MIT license (`LICENSE`, author's choice
2026-09-19); `NOTICE.md` carries the notices of the rules that descend from published ones. `Manifest.toml` files and `docs/build/` are not committed.
