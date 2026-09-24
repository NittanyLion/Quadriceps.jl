<p align="center"><img src="docs/src/assets/logo.svg" alt="Quadriceps logo" width="200"></p>

# Quadriceps.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://NittanyLion.github.io/Quadriceps.jl/dev/)
[![CI](https://github.com/NittanyLion/Quadriceps.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NittanyLion/Quadriceps.jl/actions/workflows/CI.yml)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
![authored by: JP](authored_by.svg)

> **Paper:** Joris Pinkse, *Positive weight Hermite and Legendre quadrature rules* (2026) — **[arXiv:2609.26840](https://arxiv.org/abs/2609.26840)**; Zenodo, DOI: **[10.5281/zenodo.22904159](https://doi.org/10.5281/zenodo.22904159)**
>
> **Data deposit:** Zenodo — DOI: **[10.5281/zenodo.22881864](https://doi.org/10.5281/zenodo.22881864)**

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
142 rules are stored: for example 244 nodes instead of 3125 for the Gaussian weight at
`d = 5, q = 5`, and 10984 instead of 161051 for the cube at `d = 5, q = 11`.

| `d` | GH | largest GH rule | Le | largest Le rule |
|---|---|---|---|---|
| 2 | `q ≤ 17` (`p ≤ 33`) | 208 nodes | `q ≤ 39` (`p ≤ 77`) | 1032 nodes |
| 3 | `q ≤ 17` (`p ≤ 33`) | 2226 nodes | `q ≤ 23` (`p ≤ 45`) | 4308 nodes |
| 4 | `q ≤ 12` (`p ≤ 23`) | 3238 nodes | `q ≤ 12` (`p ≤ 23`) | 3244 nodes |
| 5 | `q ≤ 11` (`p ≤ 21`) | 13199 nodes | `q ≤ 11` (`p ≤ 21`) | 10984 nodes |

<!-- END GENERATED coverage -->

The full list, with node counts, Möller's lower bound, measured accuracy and the origin of each
rule, is in the documentation's [Stored rules](https://NittanyLion.github.io/Quadriceps.jl/dev/rules.html)
page (source: [`docs/src/rules.md`](docs/src/rules.md)). These are exactly the rules of the Zenodo
deposit; the package ships nothing that is not deposited there.

## Installation

```julia
using Pkg
Pkg.add("Quadriceps")
```

Julia 1.10 or later. The only dependency is FastGaussQuadrature.jl. The package ships its rules
in double precision, in quadruple precision and to 80 digits, all in its own data files (see
below); it never reads anything from outside itself.

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
  their quadruple-precision roundings are in `data/rules128.bin` and the 80-digit rules in
  `data/rules80.bin`, with `data/index128.tsv` as the catalog of both. The
  [data format](https://NittanyLion.github.io/Quadriceps.jl/dev/format.html) page specifies all three.
* Unexported helpers: `Quadriceps.available(family)` lists the stored rules,
  `Quadriceps.nnodes(family, d, q; pragmatic)` gives a node count without building the rule,
  `Quadriceps.ruleinfo(family, d, q; pragmatic)` describes the rule and its origin (both also
  take `p` as a keyword in place of `q`), and
  `Quadriceps.exactness_error(X, w, p, family)` measures how exact a rule is. `family` is `:gh`
  or `:le`.

## Accuracy

The default data, `data/rules.bin`, are in double precision (quadruple: next section). Every
rule is its extended-precision rule rounded to `Float64`, so what remains of its error is
rounding: the largest relative monomial error over all monomials of degree `≤ p` is below
`5.2e-15` for every GH rule (most are `1e-16`–`1e-15`; the worst is `d = 5`, `p = 21`) and below
`4.5e-16` for every Le rule. All weights are positive, and all nodes of the Le rules lie strictly
inside the cube. The gate applied when the data are built, and again by the test suite on every
rule, is `1e-11`. The measured value of each rule is in the catalog (`relerr` in
`Quadriceps.ruleinfo`, and the table in `docs/src/rules.md`).

## Beyond double precision: quadruple, and 80 digits

The number type is an optional first argument of `ghpos` and `lepos`; bring the type, the
package has no dependency on it. Which data answer depends on how many significant bits the type
holds (`precision(T)`):

* **up to 113 bits** — `Float128` of Quadmath.jl, `Double64` of DoubleFloats.jl, `BigFloat` under
  `setprecision(BigFloat, 113)` — come from `data/rules128.bin`, shipped with the package: every
  rule correctly rounded to IEEE binary128 (113 bits, about 34 digits), exact for `Float128`. In
  that format the error is below `8.8e-34` for every GH rule (4.5 units in the last place; the
  worst is `d = 2`, `p = 31`) and below `3.3e-35` for every Le rule (a sixth of a unit).
* **more than 113 bits** — `BigFloat` at its default 256 bits, `Float64x4` of MultiFloats.jl, … —
  come from `data/rules80.bin`, also shipped with the package: the **80-digit rules of the Zenodo
  deposit themselves**, each number in a 40-byte binary format with a 305-bit significand. The
  result carries the deposit's 80 digits and no more, whatever the precision of the type; the
  deposit's measured error of every 80-digit rule is below `1e-68`.

All three data files are inside the package; it never reads anything from anywhere else.

```julia
using Quadriceps, Quadmath
X, w = ghpos(Float128, 3, 4)            # the 27-node rule in quadruple precision
X, w = ghpos(BigFloat, 3, 4)            # the same rule to 80 digits
X, w = lepos(BigFloat, 5; p = 21)
Float64.(X) ≈ ghpos(3, 4)[1]            # true: the same rule, rounded
```

Every stored rule is available in both ways, so the typed call never fails where the `Float64`
call succeeds; `Quadriceps.extended(family)` lists the cells with both errors (`relerr128`,
`relerr80`). Details: the
[guide](https://NittanyLion.github.io/Quadriceps.jl/dev/guide.html#Beyond-double-precision)
and the [data format](https://NittanyLion.github.io/Quadriceps.jl/dev/format.html).

## Whose rules these are

<!-- BEGIN GENERATED credits -->
116 of the 142 rules were computed from scratch by the author. A further 11 (Legendre) rules were
obtained by node elimination started from Diallo and Worku's published rules. Finally, 15 are
rules from the literature (copied in, or found again by the author's search and recognized).
<!-- END GENERATED credits -->
`Quadriceps.ruleinfo` and the
`origin` column of `data/index.tsv` say which is which; cite the source named there when you
use such a rule. [`NOTICE.md`](NOTICE.md) has the license notice that travels with the derived
files, and [`docs/src/credits.md`](docs/src/credits.md) the references.

## Documentation

The documentation is at <https://NittanyLion.github.io/Quadriceps.jl/dev/>, rebuilt by CI on
every push to `main`:

* [Guide](https://NittanyLion.github.io/Quadriceps.jl/dev/guide.html) — installation, the two
  functions, `normalize`, `pragmatic`, accuracy, number types beyond `Float64`;
* [Stored rules](https://NittanyLion.github.io/Quadriceps.jl/dev/rules.html) — every rule with its
  node count, Möller's bound, measured error and origin;
* [Data format](https://NittanyLion.github.io/Quadriceps.jl/dev/format.html) — `rules.bin`,
  `rules128.bin`, `rules80.bin` and the catalogs, enough to write a reader in any language;
* [Reference](https://NittanyLion.github.io/Quadriceps.jl/dev/api.html) — the docstrings;
* [Credits](https://NittanyLion.github.io/Quadriceps.jl/dev/credits.html) — whose rules these are.

`docs/` holds the Documenter.jl source. Build it locally with

```
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

and open `docs/build/index.html`.

## Sister packages

The same rules, with the same functions and conventions, are available for Python
([quadriceps-py](https://github.com/NittanyLion/quadriceps-py), numpy only), for R
([quadriceps-r](https://github.com/NittanyLion/quadriceps-r), base R only) and for Stata
([quadriceps-stata](https://github.com/NittanyLion/quadriceps-stata), Mata only). This package
is the master copy of the data.

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
commits and pushes, then carries the data into the three sister packages, which it tests (Stata
excepted: no Stata here), commits and pushes as well. `build/update.sh --install` adds an hourly cron entry for it (on one machine
only), `--remove` takes it out; the log is `~/.local/state/quadriceps/update.log`. Anything that
needs a person (a regression, a credit change, failing tests, a failed push) raises a desktop
notification and leaves the repositories untouched.

## License

MIT; see [`LICENSE`](LICENSE). The rules that descend from, or coincide with, published rules
carry their sources' notices in [`NOTICE.md`](NOTICE.md).
