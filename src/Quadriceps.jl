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
X, w = ghpos(3, 7)      # 27 nodes; X is 27×3, w has length 27
X, w = lepos(2, 9)      # 17 nodes on the unit square
```

Unexported but public: [`Quadriceps.available`](@ref), [`Quadriceps.nnodes`](@ref),
[`Quadriceps.ruleinfo`](@ref), [`Quadriceps.exactness_error`](@ref).
"""
module Quadriceps

using FastGaussQuadrature: gausshermite, gausslegendre

export ghpos, lepos

include("index.jl")
include("plan.jl")
include("api.jl")
include("verify.jl")

VERSION ≥ v"1.11.0-DEV.469" && eval(Meta.parse("public available, nnodes, ruleinfo, exactness_error, RuleInfo"))

end # module
