"""
    Quadriceps

Positive-weight cubature rules in `d ≥ 1` dimensions for two weight functions:

* **GH** (`ghpos`): the Gaussian weight, by default the standard normal density `N(0, I_d)`;
* **Le** (`lepos`): the uniform weight, by default the uniform density on `[0,1]^d`.

A rule of degree `p` integrates every polynomial of total degree `≤ p` exactly. All weights
are strictly positive. The rules shipped with the package are the smallest ones known to its
author; see the README and the documentation for where each one comes from.

The two exported functions mirror `gausshermite` and `gausslegendre` from
FastGaussQuadrature.jl:

```julia
X, w = ghpos(3, 4)        # as exact as the 4×4×4 Gauss–Hermite grid (degree 7), with 27 nodes
X, w = ghpos(3; p = 7)    # the same rule, requested by degree
X, w = lepos(2, 5)        # degree 9 on the unit square: 17 nodes instead of 25
```

The second argument `q` is the number of nodes of the one-dimensional Gauss rule whose
exactness is wanted, as in `gausshermite(q)`; the degree is `p = 2q - 1`.

Unexported but public: [`Quadriceps.available`](@ref), [`Quadriceps.extended`](@ref), [`Quadriceps.nnodes`](@ref),
[`Quadriceps.ruleinfo`](@ref), [`Quadriceps.exactness_error`](@ref).
"""
module Quadriceps

using FastGaussQuadrature: gausshermite, gausslegendre
using LazyArtifacts                       # the 80-digit rules, fetched from Zenodo on first use (Artifacts.toml)

export ghpos, lepos

include("index.jl")
include("extended.jl")
include("plan.jl")
include("api.jl")
include("verify.jl")

VERSION ≥ v"1.11.0-DEV.469" && eval(Meta.parse("public available, extended, nnodes, ruleinfo, exactness_error, RuleInfo"))

end # module
