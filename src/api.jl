# The exported functions.

function rule(family::Symbol, d::Integer, p::Integer, pragmatic::Bool)
    parts, n = plan(INDEX, family, d, p, pragmatic)
    materialize(family, parts, n)
end

"""
    ghpos(d, p; normalize = true, pragmatic = false) -> X, w

Positive-weight cubature rule of degree `p` for the Gaussian weight in `d` dimensions: the
smallest one the package has. `X` is an `n × d` matrix with one node per row, `w` the vector
of `n` weights, all strictly positive, and

```math
\\sum_{i=1}^n w_i f(X_{i,1}, …, X_{i,d}) = \\int f(x)\\, ω(x)\\, dx
```

for every polynomial `f` of total degree `≤ p` (up to rounding; see
[`Quadriceps.ruleinfo`](@ref) for the measured error of each stored rule).

The signature follows `gausshermite(n; normalize)` of FastGaussQuadrature.jl, with the
dimension and the degree in place of the number of nodes, and with the same meaning of
`normalize`. **The default differs**: here `normalize = true`.

# Keyword arguments

* `normalize = true`: the weight is the standard normal density
  `ω(x) = (2π)^(-d/2) exp(-|x|²/2)`, so the weights sum to 1 and the rule computes
  `E f(Z)` for `Z ~ N(0, I_d)`.
  With `normalize = false` the weight is `ω(x) = exp(-|x|²)`, the convention
  `gausshermite(n)` uses by default, and the weights sum to `π^(d/2)`.
* `pragmatic = false`: what to do when no stored rule covers `(d, p)`. With `false`, throw an
  `ArgumentError`. With `true`, return the cheapest (fewest nodes) tensor product of
  lower-dimensional rules instead: stored rules and one-dimensional Gauss–Hermite rules,
  combined over the split of `d` that minimizes the number of nodes. Such a product is a
  valid positive-weight rule of degree `p`; it is just not small. With `pragmatic = true` the
  product is also returned in the rare case that it has strictly fewer nodes than the stored
  rule, so the result is always the cheapest the package can build.

# Details

* `d = 1` returns the Gauss–Hermite rule with `p ÷ 2 + 1` nodes (as an `n × 1` matrix),
  which is optimal; it is never an error.
* Rules are stored at odd degrees. An even `p` is served by the rule for `p + 1`, and more
  generally a request is served by the smallest stored rule of degree `≥ p`.
* The result is a fresh copy; mutating it does not affect later calls. Loaded rules are
  cached, so repeated calls are cheap.

# Examples

```julia
X, w = ghpos(3, 7)                       # 27 nodes
sum(w .* X[:, 1] .^ 2 .* X[:, 2] .^ 4)   # E[Z₁² Z₂⁴] = 3

ghpos(3, 41)                             # ArgumentError: nothing stored at that degree
X, w = ghpos(3, 41; pragmatic = true)    # tensor product of lower-dimensional rules
X, w = ghpos(7, 9; pragmatic = true)     # d = 7 from a product such as (d = 4) × (d = 3)
```

See also [`lepos`](@ref), [`Quadriceps.nnodes`](@ref), [`Quadriceps.available`](@ref).
"""
function ghpos(d::Integer, p::Integer; normalize::Bool = true, pragmatic::Bool = false)
    X, w = rule(:gh, d, p, pragmatic)
    if !normalize                       # ∫ f(x) exp(-|x|²) dx = π^(d/2) E f(Z/√2)
        X ./= sqrt(2.0)
        w .*= π^(d / 2)
    end
    X, w
end

"""
    lepos(d, p; normalize = true, pragmatic = false) -> X, w

Positive-weight cubature rule of degree `p` for the uniform weight on a `d`-dimensional cube:
the smallest one the package has. `X` is an `n × d` matrix with one node per row, `w` the
vector of `n` weights, all strictly positive, and `sum(w[i] * f(X[i, :]))` equals the
integral of `f` against the weight for every polynomial `f` of total degree `≤ p` (up to
rounding; see [`Quadriceps.ruleinfo`](@ref)).

The signature follows `gausslegendre(n)` of FastGaussQuadrature.jl, with the dimension and
the degree in place of the number of nodes.

# Keyword arguments

* `normalize = true`: the weight is the uniform density on `[0,1]^d`, so the weights sum to 1
  and the rule computes `E f(U)` for `U` uniform on the unit cube.
  With `normalize = false` the rule is for `∫ f(x) dx` over `[-1,1]^d`, the convention of
  `gausslegendre`, and the weights sum to `2^d`.
* `pragmatic = false`: as for [`ghpos`](@ref). With `false`, a `(d, p)` that no stored rule
  covers is an `ArgumentError`; with `true`, the cheapest tensor product of lower-dimensional
  rules (stored rules and one-dimensional Gauss–Legendre rules) is returned instead.

The details listed under [`ghpos`](@ref) (`d = 1`, even `p`, copies and caching) apply here
too. All nodes of every stored Le rule lie inside the cube.

# Examples

```julia
X, w = lepos(2, 9)                       # 17 nodes on [0,1]²
sum(w .* X[:, 1] .^ 3 .* X[:, 2] .^ 2)   # 1/4 · 1/3

X, w = lepos(2, 9; normalize = false)    # the same rule on [-1,1]², weights sum to 4
X, w = lepos(6, 11; pragmatic = true)    # d = 6 from a product of stored rules
```

See also [`ghpos`](@ref), [`Quadriceps.nnodes`](@ref), [`Quadriceps.available`](@ref).
"""
function lepos(d::Integer, p::Integer; normalize::Bool = true, pragmatic::Bool = false)
    X, w = rule(:le, d, p, pragmatic)
    if !normalize                       # [0,1]^d → [-1,1]^d
        @. X = 2X - 1
        w .*= 2.0^d
    end
    X, w
end

"""
    Quadriceps.nnodes(family, d, p; pragmatic = false) -> Integer

Number of nodes of the rule that `ghpos(d, p; pragmatic)` (`family = :gh`) or
`lepos(d, p; pragmatic)` (`family = :le`) returns, without building it. Throws the same
`ArgumentError` when there is no rule. The count is an `Int` when it fits and a `BigInt`
otherwise (tensor products in high dimensions).
"""
function nnodes(family::Symbol, d::Integer, p::Integer; pragmatic::Bool = false)
    _, n = plan(INDEX, family, d, p, pragmatic)
    n ≤ typemax(Int) ? Int(n) : n
end

"""
    Quadriceps.ruleinfo(family, d, p; pragmatic = false)

Describe the rule behind `ghpos(d, p; pragmatic)` or `lepos(d, p; pragmatic)`: a vector with
one entry per tensor factor, in the order of the columns of `X`. An entry is either the
[`Quadriceps.RuleInfo`](@ref) of a stored rule (dimension, degree, node count, measured
error, origin) or, for a one-dimensional Gauss factor, the named tuple
`(family, d = 1, p, n, origin = "Gauss")`. A request answered by a single stored rule gives a
one-element vector.

Use it to find out whom to cite: the `origin` field names the published source of every rule
that is not the package author's own.
"""
function ruleinfo(family::Symbol, d::Integer, p::Integer; pragmatic::Bool = false)
    parts, _ = plan(INDEX, family, d, p, pragmatic)
    [a.info ≡ nothing ? (family = family, d = 1, p = 2a.n - 1, n = a.n, origin = "Gauss") : a.info for a in parts]
end
