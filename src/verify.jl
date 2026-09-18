# Exactness check of a rule against the closed-form moments of its weight.

# One-dimensional moments of the normalized weights: E Z^e for N(0,1), E U^e for U[0,1].
function moments1d(::Type{T}, family::Symbol, p::Integer) where {T}
    m = zeros(T, p + 1)
    if family ≡ :gh
        m[1] = one(T)
        for e in 2:2:p
            m[e+1] = m[e-1] * (e - 1)          # (e-1)!!
        end
    else
        for e in 0:p
            m[e+1] = one(T) / (e + 1)
        end
    end
    m
end

"""
    Quadriceps.exactness_error(X, w, p, family) -> err

Largest relative monomial error of the rule `(X, w)` over all monomials of total degree
`≤ p`, for the normalized weight of `family` (`:gh`: `N(0, I_d)`; `:le`: uniform on
`[0,1]^d`):

```math
\\max_{|a| ≤ p} \\frac{|\\sum_i w_i x_i^a - \\mathrm{E}\\,x^a|}{\\max(\\sum_i |w_i|\\,|x_i^a|,\\ 1)} .
```

The denominator is the scale of the sum being computed, so a value near `eps()` means the
rule is exact to rounding. The element type of `X` and `w` sets the arithmetic: pass
`BigFloat` arrays to measure beyond double precision. Rules returned with
`normalize = false` must be checked in the normalized frame.

```julia
X, w = ghpos(4, 5)                           # degree 2·5 - 1 = 9
Quadriceps.exactness_error(X, w, 9, :gh)     # ≈ 4e-16
```
"""
function exactness_error(X::AbstractMatrix{T}, w::AbstractVector{T}, p::Integer, family::Symbol) where {T<:AbstractFloat}
    checkfamily(family)
    n, d = size(X)
    length(w) == n || throw(DimensionMismatch("X has $n rows, w has length $(length(w))"))
    m = moments1d(T, family, p)
    P = [T[X[i, k]^e for i in 1:n, e in 0:p] for k in 1:d]      # P[k][i, e+1] = x_ik^e
    val = [Vector{T}(undef, n) for _ in 1:d]                     # val[k]: w .* x_1^a_1 ⋯ x_k^a_k
    descend(P, val, m, 1, w, one(T), Int(p))
end

# Largest error over the monomials that extend x_1^a_1 ⋯ x_(k-1)^a_(k-1) (whose weighted
# values are `prev` and whose exact moment is `mom`) by total degree ≤ r in x_k, …, x_d.
function descend(P::Vector{Matrix{T}}, val::Vector{Vector{T}}, m::Vector{T}, k::Int, prev::AbstractVector{T}, mom::T, r::Int) where {T}
    Pk = P[k]; n = length(prev); err = zero(T)
    for e in 0:r
        if k == length(P)
            s = zero(T); sa = zero(T)
            @inbounds @simd for i in 1:n
                t = prev[i] * Pk[i, e+1]
                s += t; sa += abs(t)
            end
            err = max(err, abs(s - mom * m[e+1]) / max(sa, one(T)))
        else
            v = val[k]
            @inbounds @simd for i in 1:n
                v[i] = prev[i] * Pk[i, e+1]
            end
            err = max(err, descend(P, val, m, k + 1, v, mom * m[e+1], r - e))
        end
    end
    err
end
