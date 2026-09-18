# Guide

## What is returned

```julia
X, w = ghpos(d, q)          # or ghpos(d; p = …)
X, w = lepos(d, q)          # or lepos(d; p = …)
```

`d ≥ 1` is the dimension; `q` and `p` are explained in the next section. `X` is an `n × d`
`Matrix{Float64}` with one node per row, and `w` a `Vector{Float64}` of `n` strictly positive
weights. An integral is approximated by

```julia
sum(w[i] * f(view(X, i, :)) for i in eachindex(w))
```

Every call returns fresh arrays, so the result may be modified freely. Rule files are parsed on
first use and cached, so later calls for the same rule cost a copy.

## `q` or `p`

The signatures follow `gausshermite(q; normalize)` and `gausslegendre(q)` of
FastGaussQuadrature.jl, where `q` is the number of nodes of the one-dimensional Gauss rule. That
rule is exact to degree ``2q - 1``, and so is its `d`-fold product, the ``q^d``-node grid.
`ghpos(d, q)` and `lepos(d, q)` return a rule of that same degree, ``p = 2q - 1``, with fewer
nodes than the grid. So `q` keeps its one-dimensional meaning: it says which Gauss grid the rule
stands in for. It is not the number of nodes returned, except for `d = 1`, where the result is
the `q`-node Gauss rule itself.

A second method takes the degree, as a required keyword:

```julia
ghpos(3; p = 7) == ghpos(3, 4)      # true
lepos(2; p = 12)                    # any p ≥ 0; served by the rule for p = 13
```

Rules are stored at odd degrees. A request is served by the smallest stored rule of degree
`≥ p`, since a rule exact to a higher degree is exact to degree `p`; an even `p` therefore gets
the rule for `p + 1`. The helpers [`Quadriceps.nnodes`](@ref) and
[`Quadriceps.ruleinfo`](@ref) accept `q` or `p` in the same two ways.

## Weights and the `normalize` keyword

In FastGaussQuadrature, `gausshermite(n)` integrates against ``e^{-x^2}``, and only
`gausshermite(n; normalize = true)` integrates against the standard normal density. The rules
in this package were computed for the normal density, which is what an expectation
``\mathrm{E} f(Z)`` needs, so here `normalize = true` is the default. The keyword means the
same as in FastGaussQuadrature; only the default differs.

| | `normalize = true` (default) | `normalize = false` |
|:---|:---|:---|
| [`ghpos`](@ref) | ``ω(x) = (2π)^{-d/2} e^{-\lvert x\rvert^2/2}`` on ``\mathbb{R}^d``; ``\sum w_i = 1`` | ``ω(x) = e^{-\lvert x\rvert^2}`` on ``\mathbb{R}^d``; ``\sum w_i = π^{d/2}`` |
| [`lepos`](@ref) | uniform density on ``[0,1]^d``; ``\sum w_i = 1`` | ``ω = 1`` on ``[-1,1]^d``; ``\sum w_i = 2^d`` |

With `normalize = false` both functions follow FastGaussQuadrature's default conventions, so
that `ghpos(1, q; normalize = false)` is `gausshermite(q)` and
`lepos(1, q; normalize = false)` is `gausslegendre(q)`, up to the shape of `X`. The two
frames are related by

```math
\text{GH:}\quad x ↦ x/\sqrt{2},\;\; w ↦ π^{d/2} w ;
\qquad
\text{Le:}\quad x ↦ 2x - 1,\;\; w ↦ 2^d w .
```

Other Gaussian and rectangular weights follow by a change of variables. For
``Y \sim N(μ, Σ)`` with ``Σ = LL^\top``, the nodes ``μ + L x_i`` with the same weights give a
rule of the same degree for ``Y``:

```julia
X, w = ghpos(3, 5)
Y = μ' .+ X * L'            # rows are the nodes for N(μ, LL')
```

For the uniform density on a box ``\prod_k [a_k, b_k]``, use the nodes ``a + (b - a) ⊙ x_i``
with the same weights.

## Requests without a stored rule: the `pragmatic` keyword

Rules are stored for ``2 ≤ d ≤ 5`` and odd degrees up to a ceiling that depends on the family
and on ``d`` (see [Stored rules](rules.md)). One dimension needs no storage: the Gauss rule
is optimal and is computed on the fly.

When no stored rule of dimension `d` has degree `≥ p`:

* with `pragmatic = false`, the default, the functions throw an `ArgumentError` whose message
  says how far the stored rules for that dimension go;
* with `pragmatic = true` they return the **cheapest tensor product of lower-dimensional
  rules**.

A tensor product of rules of degree `p` for the factors of a product weight is a rule of degree
`p` for the product weight, with positive weights if the factors have them; both weights here
are product weights. The package considers every way of writing ``d = d_1 + ⋯ + d_k`` with
each factor either a stored rule of dimension ``d_j`` and degree `≥ p` or a one-dimensional
Gauss rule, and returns the combination with the fewest nodes. The columns of `X` follow the
order of the factors, and the first factor varies slowest along the rows.

```@repl guide
using Quadriceps
ghpos(7, 5)
X, w = ghpos(7, 5; pragmatic = true); size(X)
[(r.d, r.n) for r in Quadriceps.ruleinfo(:gh, 7, 5; pragmatic = true)]
5^7     # the product grid
```

When the degree, not the dimension, is out of range, the saving is smaller, because a
one-dimensional Gauss factor is unavoidable:

```@repl guide
[(r.d, r.n) for r in Quadriceps.ruleinfo(:le, 3, 24; pragmatic = true)]
Quadriceps.nnodes(:le, 3, 24; pragmatic = true), 24^3
```

`pragmatic = true` changes nothing for a request that a stored rule covers, with one proviso:
the result is always the cheapest rule the package can build, so if a tensor product had
strictly fewer nodes than the stored rule, the product would be returned. The data contain no
such case. When the data are built, a rule that a tensor product matches is left out, unless
it attains Möller's bound; the test suite checks this.

Node counts grow quickly. [`Quadriceps.nnodes`](@ref) returns the count without building the
rule, and a request whose rule would not fit in memory is refused with an `ArgumentError`.

## Accuracy

The rules are stored in double precision, with 17 significant digits. When the package data
are built, every rule is checked:

* all weights are strictly positive;
* the largest relative monomial error over all monomials of total degree `≤ p`,
  ```math
  \max_{|a| ≤ p} \frac{|\sum_i w_i x_i^a - \mathrm{E}\,x^a|}{\max(\sum_i |w_i|\,|x_i^a|,\ 1)} ,
  ```
  is below ``10^{-11}``.

That threshold is a gate, not the accuracy achieved. The measured error of each rule is
recorded in the catalog and shown in [Stored rules](rules.md): most are between ``10^{-16}``
and ``10^{-15}``. The largest GH rules for ``d = 2`` (``p ≥ 33``) and ``d = 3`` (``p ≥ 31``)
are the least accurate, between ``10^{-12}`` and ``10^{-11}``; their error is a property of the
stored nodes and weights, not of the arithmetic of the check.
[`Quadriceps.exactness_error`](@ref) repeats the measurement, in extended precision if given
`BigFloat` arrays, and the test suite runs it on every stored rule.

All nodes of the Le rules lie strictly inside the cube. GH nodes are unrestricted.

A rule with `n` equal to Möller's lower bound is proven minimal: no rule of that degree, with
or without positive weights, has fewer nodes. Those are marked in [Stored rules](rules.md).
Everywhere else "smallest" means smallest known to the author, not smallest possible.

## Inspecting the catalog

```@repl guide
using Quadriceps: available, nnodes, ruleinfo
[(r.p, r.n) for r in available(:le) if r.d == 4]     # degrees and node counts, Le, d = 4
nnodes(:gh, 5, 7)                                     # q = 7, degree 13
ruleinfo(:le, 3; p = 41)[1].origin
```

All rules are stored in one binary file, `data/rules.bin`; `data/index.tsv` is the catalog that
goes with it. [Data format](format.md) specifies both.
