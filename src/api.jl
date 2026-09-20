# The exported functions.

# Degree of the q-node Gauss rule, which is what a request for q asks to match.
function degree(q::Integer)
    q ≥ 1 || throw(ArgumentError("q must be at least 1, got $q"))
    2q - 1
end

function rule(::Type{T}, family::Symbol, d::Integer, p::Integer, pragmatic::Bool) where {T<:AbstractFloat}
    parts, n = plan(INDEX, family, d, p, pragmatic)
    materialize(T, family, parts, n)
end

"""
    ghpos([T = Float64,] d, q; normalize = true, pragmatic = false) -> X, w
    ghpos([T = Float64,] d; p, normalize = true, pragmatic = false) -> X, w

Positive-weight cubature rule for the Gaussian weight in `d` dimensions: the smallest one the
package has.

The first method follows `gausshermite(q; normalize)` of FastGaussQuadrature.jl. There, `q` is
the number of nodes of the one-dimensional Gauss rule, which is exact to degree `2q - 1`. Here,
`ghpos(d, q)` returns a rule with the same exactness in `d` dimensions, degree `p = 2q - 1`,
that replaces the `q^d`-node product grid; for `d = 1` it is `gausshermite(q)` itself. The
second method takes the degree directly, as the keyword `p`; `ghpos(d, q)` is
`ghpos(d; p = 2q - 1)`.

`X` is an `n × d` matrix with one node per row, `w` the vector of `n` weights, all strictly
positive, and

```math
\\sum_{i=1}^n w_i f(X_{i,1}, …, X_{i,d}) = \\int f(x)\\, ω(x)\\, dx
```

for every polynomial `f` of total degree `≤ p` (up to rounding; see
[`Quadriceps.ruleinfo`](@ref) for the measured error of each stored rule).

# Keyword arguments

* `p`: the degree of exactness, `p ≥ 0` (second method only, required). Rules are stored at
  odd degrees; a request is served by the smallest stored rule of degree `≥ p`, so an even `p`
  gets the rule for `p + 1`.
* `normalize = true`: the weight is the standard normal density
  `ω(x) = (2π)^(-d/2) exp(-|x|²/2)`, so the weights sum to 1 and the rule computes
  `E f(Z)` for `Z ~ N(0, I_d)`.
  With `normalize = false` the weight is `ω(x) = exp(-|x|²)` and the weights sum to
  `π^(d/2)`. The keyword means what it means in `gausshermite`, but **the default differs**:
  `gausshermite` defaults to `normalize = false`, and only its `normalize = true` is a rule for
  the normal density.
* `pragmatic = false`: what to do when no stored rule covers the request. With `false`, throw
  an `ArgumentError`. With `true`, return the cheapest (fewest nodes) tensor product of
  lower-dimensional rules instead: stored rules and one-dimensional Gauss–Hermite rules,
  combined over the split of `d` that minimizes the number of nodes. Such a product is a
  valid positive-weight rule of the requested degree; it is just not small. With
  `pragmatic = true` the product is also returned in the rare case that it has strictly fewer
  nodes than the stored rule, so the result is always the cheapest the package can build.

# Number type

The optional first argument `T` is the element type of `X` and `w`. Beyond `Float64` — say
`Float128` of Quadmath.jl, `Double64` of DoubleFloats.jl, or `BigFloat` — the rule comes from
the package's quadruple-precision data: every node and weight is the stored extended-precision
rule correctly rounded to IEEE binary128 (113 bits, about 34 digits); rounded on to `Float64` it
is the `Float64` rule (to the last bit, except that a coordinate of size `1e-30` standing for
zero can come out one unit in the last place off). A `BigFloat` result therefore carries 34 correct
digits, not more. One-dimensional Gauss factors are computed to that accuracy as well.
[`Quadriceps.extended`](@ref) lists the cells stored this way, with the measured error of each
in that format: every stored rule is one of them, so `T` never fails on a cell the `Float64` call
answers.

```julia
using Quadmath
X, w = ghpos(Float128, 3, 4)             # the 27-node rule in quadruple precision
Float64.(X) ≈ ghpos(3, 4)[1]             # true: the same rule
```

# Details

* `d = 1` returns the Gauss–Hermite rule with `q` (or `p ÷ 2 + 1`) nodes, as an `n × 1`
  matrix; it is never an error.
* The result is a fresh copy; mutating it does not affect later calls. Loaded rules are
  cached, so repeated calls are cheap.

# Examples

```julia
X, w = ghpos(3, 4)                       # degree 7: 27 nodes instead of 4³ = 64
sum(w .* X[:, 1] .^ 2 .* X[:, 2] .^ 4)   # E[Z₁² Z₂⁴] = 3
ghpos(3; p = 7) == ghpos(3, 4)           # true

ghpos(3, 21)                             # ArgumentError: nothing stored at that degree
X, w = ghpos(3, 21; pragmatic = true)    # tensor product of lower-dimensional rules
X, w = ghpos(7, 5; pragmatic = true)     # d = 7 as (d = 2) × (d = 5)
```

See also [`lepos`](@ref), [`Quadriceps.nnodes`](@ref), [`Quadriceps.available`](@ref).
"""
ghpos(d::Integer, q::Integer; kw...) = ghpos(Float64, d, q; kw...)
ghpos(d::Integer; kw...) = ghpos(Float64, d; kw...)
ghpos(::Type{T}, d::Integer, q::Integer; normalize::Bool = true, pragmatic::Bool = false) where {T<:AbstractFloat} =
    ghpos(T, d; p = degree(q), normalize, pragmatic)

function ghpos(::Type{T}, d::Integer; p::Integer, normalize::Bool = true, pragmatic::Bool = false) where {T<:AbstractFloat}
    X, w = rule(T, :gh, d, p, pragmatic)
    if !normalize                       # ∫ f(x) exp(-|x|²) dx = π^(d/2) E f(Z/√2)
        X ./= sqrt(T(2))
        w .*= T(π)^(T(d) / 2)
    end
    X, w
end

"""
    lepos([T = Float64,] d, q; normalize = true, pragmatic = false) -> X, w
    lepos([T = Float64,] d; p, normalize = true, pragmatic = false) -> X, w

Positive-weight cubature rule for the uniform weight on a `d`-dimensional cube: the smallest
one the package has.

The first method follows `gausslegendre(q)` of FastGaussQuadrature.jl: `q` is the number of
nodes of the one-dimensional Gauss rule, and `lepos(d, q)` returns a rule with the same
exactness in `d` dimensions, degree `p = 2q - 1`, that replaces the `q^d`-node product grid;
for `d = 1` it is `gausslegendre(q)` itself. The second method takes the degree directly, as
the keyword `p`; `lepos(d, q)` is `lepos(d; p = 2q - 1)`.

`X` is an `n × d` matrix with one node per row, `w` the vector of `n` weights, all strictly
positive, and `sum(w[i] * f(X[i, :]))` equals the integral of `f` against the weight for every
polynomial `f` of total degree `≤ p` (up to rounding; see [`Quadriceps.ruleinfo`](@ref)).

# Keyword arguments

* `p`: the degree of exactness, `p ≥ 0` (second method only, required); an even `p` gets the
  rule for `p + 1`.
* `normalize = true`: the weight is the uniform density on `[0,1]^d`, so the weights sum to 1
  and the rule computes `E f(U)` for `U` uniform on the unit cube.
  With `normalize = false` the rule is for `∫ f(x) dx` over `[-1,1]^d`, the convention of
  `gausslegendre`, and the weights sum to `2^d`.
* `pragmatic = false`: as for [`ghpos`](@ref). With `false`, a request that no stored rule
  covers is an `ArgumentError`; with `true`, the cheapest tensor product of lower-dimensional
  rules (stored rules and one-dimensional Gauss–Legendre rules) is returned instead.

The optional number type `T` and the details listed under [`ghpos`](@ref) apply here too. All nodes of every stored Le rule lie
inside the cube.

# Examples

```julia
X, w = lepos(2, 5)                       # degree 9: 17 nodes on [0,1]² instead of 25
sum(w .* X[:, 1] .^ 3 .* X[:, 2] .^ 2)   # 1/4 · 1/3

X, w = lepos(2, 5; normalize = false)    # the same rule on [-1,1]², weights sum to 4
X, w = lepos(6, 6; pragmatic = true)     # d = 6 from a product of stored rules
```

See also [`ghpos`](@ref), [`Quadriceps.nnodes`](@ref), [`Quadriceps.available`](@ref).
"""
lepos(d::Integer, q::Integer; kw...) = lepos(Float64, d, q; kw...)
lepos(d::Integer; kw...) = lepos(Float64, d; kw...)
lepos(::Type{T}, d::Integer, q::Integer; normalize::Bool = true, pragmatic::Bool = false) where {T<:AbstractFloat} =
    lepos(T, d; p = degree(q), normalize, pragmatic)

function lepos(::Type{T}, d::Integer; p::Integer, normalize::Bool = true, pragmatic::Bool = false) where {T<:AbstractFloat}
    X, w = rule(T, :le, d, p, pragmatic)
    if !normalize                       # [0,1]^d → [-1,1]^d
        @. X = 2X - 1
        w .*= T(2)^d
    end
    X, w
end

"""
    Quadriceps.nnodes(family, d, q; pragmatic = false) -> Integer
    Quadriceps.nnodes(family, d; p, pragmatic = false) -> Integer

Number of nodes of the rule that [`ghpos`](@ref) (`family = :gh`) or [`lepos`](@ref)
(`family = :le`) returns for the same arguments, without building it. Throws the same
`ArgumentError` when there is no rule. The count is an `Int` when it fits and a `BigInt`
otherwise (tensor products in high dimensions).
"""
nnodes(family::Symbol, d::Integer, q::Integer; pragmatic::Bool = false) = nnodes(family, d; p = degree(q), pragmatic)

function nnodes(family::Symbol, d::Integer; p::Integer, pragmatic::Bool = false)
    _, n = plan(INDEX, family, d, p, pragmatic)
    n ≤ typemax(Int) ? Int(n) : n
end

"""
    Quadriceps.ruleinfo(family, d, q; pragmatic = false)
    Quadriceps.ruleinfo(family, d; p, pragmatic = false)

Describe the rule that [`ghpos`](@ref) or [`lepos`](@ref) returns for the same arguments: a
vector with one entry per tensor factor, in the order of the columns of `X`. An entry is
either the [`Quadriceps.RuleInfo`](@ref) of a stored rule (dimension, degree, node count,
measured error, origin) or, for a one-dimensional Gauss factor, the named tuple
`(family, d = 1, p, n, origin = "Gauss")`. A request answered by a single stored rule gives a
one-element vector.

Use it to find out whom to cite: the `origin` field names the published source of every rule
that is not the package author's own.
"""
ruleinfo(family::Symbol, d::Integer, q::Integer; pragmatic::Bool = false) = ruleinfo(family, d; p = degree(q), pragmatic)

function ruleinfo(family::Symbol, d::Integer; p::Integer, pragmatic::Bool = false)
    parts, _ = plan(INDEX, family, d, p, pragmatic)
    [a.info ≡ nothing ? (family = family, d = 1, p = 2a.n - 1, n = a.n, origin = "Gauss") : a.info for a in parts]
end
