# Quadriceps.jl

Positive-weight cubature rules in several dimensions, for two weights:

| function | weight (default) | one-dimensional cousin |
|---|---|---|
| `ghpos(d, p)` | standard normal density `N(0, I_d)` on `R^d` | `gausshermite` |
| `lepos(d, p)` | uniform density on `[0,1]^d` | `gausslegendre` |

A rule of degree `p` is a set of `n` nodes `x_i ∈ R^d` and weights `w_i > 0` with
`Σ w_i f(x_i) = ∫ f(x) ω(x) dx` for every polynomial `f` of total degree `≤ p`. The product of
one-dimensional Gauss rules does this with `q^d` nodes, `q = (p+1)/2`. The rules stored here do
it with far fewer: 244 nodes instead of 3125 for the Gaussian weight at `d = 5, p = 9`, 10984
instead of 161051 for the cube at `d = 5, p = 21`. They are the smallest positive-weight rules
known to the author, 148 in all:

| `d` | GH degrees | largest GH rule | Le degrees | largest Le rule |
|---|---|---|---|---|
| 2 | 1, 3, …, 43 | 482 nodes | 1, 3, …, 77 | 1032 nodes |
| 3 | 1, 3, …, 35 | 4750 nodes | 1, 3, …, 45 | 4308 nodes |
| 4 | 1, 3, …, 23 | 3238 nodes | 1, 3, …, 23 | 3244 nodes |
| 5 | 1, 3, …, 21 | 13199 nodes | 1, 3, …, 21 | 10984 nodes |

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

X, w = ghpos(3, 7)              # d = 3, degree 7: X is 27×3 (one node per row), w has length 27
f(x) = x[1]^2 * x[2]^4
sum(w[i] * f(X[i, :]) for i in eachindex(w))        # E[Z₁² Z₂⁴] = 3.0

X, w = lepos(2, 9)              # 17 nodes on the unit square
sum(w .* X[:, 1] .^ 3 .* X[:, 2] .^ 2)              # 1/4 · 1/3
```

Both functions take the dimension `d ≥ 1` and the degree `p ≥ 0` and return `(X, w)`, in the
manner of `gausshermite(n)` and `gausslegendre(n)` from FastGaussQuadrature.jl, with `(d, p)` in
place of the number of nodes. `X` is an `n × d` matrix, `w` a vector of `n` positive weights.

### `normalize`

FastGaussQuadrature's `gausshermite(n)` integrates against `exp(-x²)`; only
`gausshermite(n; normalize = true)` gives a rule for the standard normal density. The rules here
are made for the normal density, so **`normalize = true` is the default** — the opposite of
FastGaussQuadrature's default, with the same meaning of the keyword:

| | `normalize = true` (default) | `normalize = false` (FastGaussQuadrature's convention) |
|---|---|---|
| `ghpos` | weight `(2π)^(-d/2) exp(-|x|²/2)`; weights sum to 1 | weight `exp(-|x|²)`; weights sum to `π^(d/2)` |
| `lepos` | uniform density on `[0,1]^d`; weights sum to 1 | `∫ f(x) dx` over `[-1,1]^d`; weights sum to `2^d` |

So `ghpos(1, 2n-1; normalize = false)` is `gausshermite(n)` and `lepos(1, 2n-1; normalize = false)`
is `gausslegendre(n)`, up to the shape of `X`.

### `pragmatic`

Rules are stored for `2 ≤ d ≤ 5` up to the degrees in the table above. For any other `(d, p)`:

* `pragmatic = false` (the default) throws an `ArgumentError`, which says how far the stored
  rules go;
* `pragmatic = true` returns the cheapest tensor product of lower-dimensional rules: the split
  of `d` into stored rules and one-dimensional Gauss rules that needs the fewest nodes. The
  result is a valid positive-weight rule of degree `p`. It is not small, but it is much
  smaller than the plain product grid whenever a stored rule can be a factor.

```julia
ghpos(7, 9)                             # ArgumentError: no stored rule in seven dimensions
X, w = ghpos(7, 9; pragmatic = true)    # 4392 nodes: (d = 2, n = 18) × (d = 5, n = 244); the grid has 78125
X, w = lepos(3, 47; pragmatic = true)   # 9312 nodes: Gauss (24) × (d = 2, n = 388); the grid has 13824

Quadriceps.nnodes(:gh, 10, 5; pragmatic = true)     # 1024, without building the rule
Quadriceps.ruleinfo(:gh, 7, 9; pragmatic = true)    # the factors, with their origins
```

With `pragmatic = true` a request that a stored rule covers returns that stored rule, as without
it. (The one exception: if a product were ever strictly cheaper than the stored rule, the
product is returned. No such case exists in the current data; the build leaves out any rule
that a product matches.)

### Other details

* `d = 1` gives the Gauss rule with `p ÷ 2 + 1` nodes, as an `n × 1` matrix.
* Rules are stored at odd degrees. A request is served by the smallest stored rule of degree
  `≥ p`, so an even `p` gets the rule for `p + 1`.
* Every call returns fresh arrays. Rule files are parsed on first use and cached.
* Unexported helpers: `Quadriceps.available(family)` lists the stored rules,
  `Quadriceps.nnodes(family, d, p; pragmatic)` gives a node count without building the rule,
  `Quadriceps.ruleinfo(family, d, p; pragmatic)` describes the rule and its origin, and
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

122 of the 148 rules were computed by the author. 15 are rules from the literature (copied in,
or found again by the author's search and recognized), and 11 Le rules were obtained by node
elimination started from Diallo and Worku's published rules. `Quadriceps.ruleinfo` and the
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

## Rebuilding the data

`build/build_data.jl` regenerates `data/` and `docs/src/rules.md` from the rule bank of the
designed-quadrature project: for every `(d, p)` it takes the smallest rule that passes the
checks, and keeps it if it has fewer nodes than every tensor product of lower-dimensional rules
or attains Möller's bound.

```
julia --project=build build/build_data.jl [path to the project's sync folder]
```
