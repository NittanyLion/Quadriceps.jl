# Quadriceps.jl

[![CI](https://github.com/NittanyLion/Quadriceps.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NittanyLion/Quadriceps.jl/actions/workflows/CI.yml)

> **Paper:** J. Pinkse, *Positive weight Hermite and Legendre quadrature rules* — arXiv: **[ARXIV-LINK-TBA](https://arxiv.org/abs/ARXIV-LINK-TBA)** (link to be filled in on publication)
>
> **Data deposit:** Zenodo — DOI: **[ZENODO-DOI-TBA](https://doi.org/ZENODO-DOI-TBA)** (link to be filled in on publication)

Positive-weight cubature rules in several dimensions, for two weights:

| function | weight (default) | one-dimensional cousin |
|---|---|---|
| `ghpos(d, q)` | standard normal density `N(0, I_d)` on `ℝᵈ` | `gausshermite(q)` |
| `lepos(d, q)` | uniform density on `[0,1]ᵈ` | `gausslegendre(q)` |

A rule of degree `p` is a set of `n` nodes `x_i ∈ ℝᵈ` and weights `w_i > 0` with
`Σ w_i f(x_i) = ∫ f(x) ω(x) dx` for every polynomial `f` of total degree `≤ p`. The product of
`q`-node one-dimensional Gauss rules does this for `p = 2q - 1` with `qᵈ` nodes. The rules
stored here, the smallest positive-weight rules known to the author, do it with far fewer.

<!-- BEGIN GENERATED coverage -->
148 rules are stored: for example 244 nodes instead of 3125 for the Gaussian weight at
`d = 5, q = 5`, and 10984 instead of 161051 for the cube at `d = 5, q = 11`.

| `d` | GH | largest GH rule | Le | largest Le rule |
|---|---|---|---|---|
| 2 | `q ≤ 22` (`p ≤ 43`) | 482 nodes | `q ≤ 39` (`p ≤ 77`) | 1032 nodes |
| 3 | `q ≤ 18` (`p ≤ 35`) | 4749 nodes | `q ≤ 23` (`p ≤ 45`) | 4308 nodes |
| 4 | `q ≤ 12` (`p ≤ 23`) | 3238 nodes | `q ≤ 12` (`p ≤ 23`) | 3244 nodes |
| 5 | `q ≤ 11` (`p ≤ 21`) | 13199 nodes | `q ≤ 11` (`p ≤ 21`) | 10984 nodes |

<!-- END GENERATED coverage -->

The full list, with node counts, Möller's lower bound, measured accuracy and the origin of each
rule, is in [`docs/src/rules.md`](docs/src/rules.md).

## Installation

The repository is private. With access to it:

```julia
using Pkg
Pkg.add(url = "git@github.com:NittanyLion/Quadriceps.jl.git")
```

Julia 1.10 or later. The only dependency is FastGaussQuadrature.jl.

## Use

```julia
using Quadriceps

X, w = ghpos(3, 4)              # d = 3, q = 4 (degree 7): X is 27×3 (one node per row), w has length 27
f(x) = x[1]^2 * x[2]^4
sum(w[i] * f(X[i, :]) for i in eachindex(w))        # E[Z₁² Z₂⁴] = 3.0

X, w = lepos(2, 5)              # q = 5 (degree 9): 17 nodes on the unit square instead of 25
sum(w .* X[:, 1] .^ 3 .* X[:, 2] .^ 2)              # 1/4 · 1/3

ghpos(3; p = 7) == ghpos(3, 4)  # true: the keyword p requests a rule by its degree
```

Both functions follow `gausshermite(q)` and `gausslegendre(q)` from FastGaussQuadrature.jl, with
the dimension `d ≥ 1` in front. There, `q` is the number of nodes of the one-dimensional Gauss
rule, which is exact to degree `2q - 1`. Here, `ghpos(d, q)` returns a `d`-dimensional rule of
that same degree `p = 2q - 1`: a replacement for the `qᵈ`-node product grid, and for `d = 1`
the `q`-node Gauss rule itself. `X` is an `n × d` matrix, `w` a vector of `n` positive weights.

To ask for a degree instead, use the method with the keyword `p`: `ghpos(d; p = 7)`,
`lepos(d; p = 12)`. Any `p ≥ 0` is accepted. Rules are stored at odd degrees and a request is
served by the smallest stored rule of degree `≥ p`, so an even `p` gets the rule for `p + 1`.

### `normalize`

FastGaussQuadrature's `gausshermite(n)` integrates against `exp(-x²)`; only
`gausshermite(n; normalize = true)` gives a rule for the standard normal density. The rules here
are made for the normal density, so **`normalize = true` is the default** — the opposite of
FastGaussQuadrature's default, with the same meaning of the keyword:

| | `normalize = true` (default) | `normalize = false` (FastGaussQuadrature's convention) |
|---|---|---|
| `ghpos` | weight `(2π)⁻ᵈᐟ² exp(-‖x‖²/2)`; weights sum to 1 | weight `exp(-‖x‖²)`; weights sum to `πᵈᐟ²` |
| `lepos` | uniform density on `[0,1]ᵈ`; weights sum to 1 | `∫ f(x) dx` over `[-1,1]ᵈ`; weights sum to `2ᵈ` |

So `ghpos(1, q; normalize = false)` is `gausshermite(q)`, `ghpos(1, q)` is
`gausshermite(q; normalize = true)`, and `lepos(1, q; normalize = false)` is `gausslegendre(q)`,
up to the shape of `X`.

### `pragmatic`

Rules are stored for `2 ≤ d ≤ 5` up to the ceilings in the table above. For any other request:

* `pragmatic = false` (the default) throws an `ArgumentError`, which says how far the stored
  rules go;
* `pragmatic = true` returns the cheapest tensor product of lower-dimensional rules: the split
  of `d` into stored rules and one-dimensional Gauss rules that needs the fewest nodes. The
  result is a valid positive-weight rule of degree `p`. It is not small, but it is much
  smaller than the plain product grid whenever a stored rule can be a factor.

<!-- BEGIN GENERATED pragmatic -->
```julia
ghpos(7, 5)                             # ArgumentError: no stored rule in seven dimensions
X, w = ghpos(7, 5; pragmatic = true)    # 4392 nodes; the product grid has 78125
X, w = lepos(3, 24; pragmatic = true)   # 9312 nodes; the product grid has 13824

Quadriceps.nnodes(:gh, 10, 3; pragmatic = true)     # 1024, without building the rule
Quadriceps.ruleinfo(:gh, 7, 5; pragmatic = true)    # the factors, with their origins
```
<!-- END GENERATED pragmatic -->

With `pragmatic = true` a request that a stored rule covers returns that stored rule, as without
it. (The one exception: if a product were ever strictly cheaper than the stored rule, the
product is returned. No such case exists in the current data; the build leaves out any rule
that a product matches.)

### Other details

* `d = 1` gives the Gauss rule with `q` nodes (`p ÷ 2 + 1` when `p` is given), as an `n × 1`
  matrix.
* Every call returns fresh arrays. A rule is read from the data file on first use and cached.
* All rules live in one binary file, `data/rules.bin`, with `data/index.tsv` as its catalog;
  [`docs/src/format.md`](docs/src/format.md) specifies the format.
* Unexported helpers: `Quadriceps.available(family)` lists the stored rules,
  `Quadriceps.nnodes(family, d, q; pragmatic)` gives a node count without building the rule,
  `Quadriceps.ruleinfo(family, d, q; pragmatic)` describes the rule and its origin (both also
  take `p` as a keyword in place of `q`), and
  `Quadriceps.exactness_error(X, w, p, family)` measures how exact a rule is. `family` is `:gh`
  or `:le`.

## Accuracy

Rules are stored in double precision. Every stored rule was checked when the data were built:
all weights positive, and the largest relative monomial error over all monomials of degree
`≤ p` below `1e-11`. Most rules sit at `1e-16`–`1e-15`; the largest GH rules at `d = 2, 3` are
the least accurate, at `1e-12`–`1e-11`. The measured value of each rule is in the catalog
(`relerr` in `Quadriceps.ruleinfo`, and the table in `docs/src/rules.md`), and the test suite
repeats the check for every rule. All nodes of the Le rules lie strictly inside the cube.

## Whose rules these are

<!-- BEGIN GENERATED credits -->
122 of the 148 rules were computed by the author. 15 are rules from the literature (copied in,
or found again by the author's search and recognized), and 11 Le rules were obtained by node
elimination started from Diallo and Worku's published rules.
<!-- END GENERATED credits -->
`Quadriceps.ruleinfo` and the
`origin` column of `data/index.tsv` say which is which; cite the source named there when you
use such a rule. [`NOTICE.md`](NOTICE.md) has the license notice that travels with the derived
files, and [`docs/src/credits.md`](docs/src/credits.md) the references.

## Documentation

`docs/` holds a Documenter.jl site: a guide, the API reference, the table of stored rules and
the credits. Build it locally with

```
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

and open `docs/build/index.html`.

## Sister packages

The same rules, with the same functions and conventions, are available for Python
([quadriceps-py](https://github.com/NittanyLion/quadriceps-py), numpy only) and for R
([quadriceps-r](https://github.com/NittanyLion/quadriceps-r), base R only). This package is the
master copy of the data.

## Keeping the data current

`build/build_data.jl` regenerates `data/`, `docs/src/rules.md` and the generated blocks of this
README from the rule bank of the designed-quadrature project: for every `(d, p)` it takes the
smallest rule that passes the checks, and keeps it if it has fewer nodes than every tensor
product of lower-dimensional rules or attains Möller's bound. It refuses to replace `data/` if a
stored rule would disappear or grow, or if the set of third-party rules would change.

```
julia --project=build build/build_data.jl [--force] [path to the project's sync folder]
```

`build/update.sh` does this unattended: it rebuilds, and when a rule changed it runs the tests,
commits and pushes, then carries the data into the two sister packages, which it tests, commits
and pushes as well. `build/update.sh --install` adds an hourly cron entry for it (on one machine
only), `--remove` takes it out; the log is `~/.local/state/quadriceps/update.log`. Anything that
needs a person (a regression, a credit change, failing tests, a failed push) raises a desktop
notification and leaves the repositories untouched.
