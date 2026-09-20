# Rules beyond double precision: data/rules128.bin and the number-type argument of ghpos / lepos.
#
# data/rules128.bin is a QUADRICEPS1 file with float=binary128: every number is the IEEE
# binary128 (quadruple precision, 113-bit significand) rounding of the project's
# extended-precision rule, 16 little-endian bytes. It is decoded here into 113-bit BigFloats, so
# the package needs no quadruple-precision dependency; `T.(…)` then gives the caller's type, and
# for `Float128` of Quadmath.jl that conversion is exact. data/index128.tsv is its catalog.

const BIN128 = BinIndex()
const INDEX128 = Dict{Tuple{Symbol,Int,Int},Float64}()          # key => relerr128
const CACHE128 = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{BigFloat},Vector{BigFloat}}}()
const BIAS128 = 16383

bin128path() = joinpath(datadir(), "rules128.bin")

# The binary128 bit pattern of a BigFloat that holds at most 113 significant bits.
function encode128(x::BigFloat)
    iszero(x) && return UInt128(0)
    e = exponent(x)
    -16382 ≤ e ≤ 16383 || throw(ArgumentError("$x is outside the normal range of binary128"))
    m = ldexp(abs(x), 112 - e)                                   # in [2^112, 2^113)
    isinteger(m) || throw(ArgumentError("$x has more than 113 significant bits"))
    (UInt128(signbit(x)) << 127) | (UInt128(e + BIAS128) << 112) | (UInt128(BigInt(m)) - (UInt128(1) << 112))
end

function decode128(b::UInt128)
    iszero(b << 1) && return BigFloat(0; precision = 113)
    e = Int((b >> 112) & 0x7fff)
    0 < e < 0x7fff || error("rules128.bin holds a subnormal or non-finite number")
    m = (b & ((UInt128(1) << 112) - 1)) | (UInt128(1) << 112)
    x = ldexp(BigFloat(BigInt(m); precision = 113), e - BIAS128 - 112)
    iszero(b >> 127) ? x : -x
end

function readblock128(path::AbstractString, e::BinEntry, d::Integer)
    v = Vector{UInt128}(undef, e.n * (d + 1))
    open(path) do io
        seek(io, e.offset)
        read!(io, v)
    end
    A = permutedims(reshape(decode128.(ltoh.(v)), d + 1, e.n))   # row-major on disk
    A[:, 1:d], A[:, d+1]
end

function readindex128(path::AbstractString)
    out = Dict{Tuple{Symbol,Int,Int},Float64}()
    isfile(path) || return out
    for l in eachline(path)
        (isempty(l) || startswith(l, '#') || startswith(l, "family")) && continue
        r = split(l, '\t')
        out[(Symbol(r[1]), parse(Int, r[2]), parse(Int, r[3]))] = parse(Float64, r[5])
    end
    out
end

# The stored rule in extended precision (113-bit BigFloats, shared: callers convert, which copies).
# Since the 2026-09-19 QUAD_ONLY policy the catalog holds no cell without a binary128 rule (the
# builder leaves such a cell out), so the error below is a guard against damaged data, not a case
# a caller can reach through ghpos / lepos.
function stored128(info::RuleInfo)
    key = (info.family, info.d, info.p)
    (haskey(BIN128, key) && BIN128[key].n == info.n) ||
        throw(ArgumentError("the rule for family :$(info.family), d = $(info.d), p = $(info.p) is stored in Float64 only; " *
                            "see Quadriceps.extended(:$(info.family)) for the cells held beyond double precision"))
    lock(CACHE_LOCK) do
        get!(CACHE128, key) do
            readblock128(bin128path(), BIN128[key], info.d)
        end
    end
end

# A number type that Float64 data serve in full.
narrow(::Type{T}) where {T} = T ≡ Float64 || T ≡ Float32 || T ≡ Float16

function stored(::Type{T}, info::RuleInfo) where {T<:AbstractFloat}
    X, w = narrow(T) ? stored(info) : stored128(info)
    T.(X), T.(w)
end

# One-dimensional Gauss rules in T: the Float64 nodes of FastGaussQuadrature.jl, refined by
# Newton's method on the orthonormal polynomial in 256-bit arithmetic; weights from the
# Christoffel function, 1 / Σ_{k<q} φ_k(x)², which sum to 1 for a probability weight.
function gauss1d(::Type{T}, family::Symbol, q::Integer) where {T<:AbstractFloat}
    narrow(T) && return (r = gauss1d(family, q); (T.(r[1]), T.(r[2])))
    x, w = setprecision(BigFloat, 256) do
        x = BigFloat.(family ≡ :gh ? gausshermite(q; normalize = true)[1] : gausslegendre(q)[1])
        w = similar(x)
        for i in eachindex(x)
            for _ in 1:6
                φ, dφ, _ = orthonormal(family, q, x[i])
                x[i] -= φ / dφ
            end
            w[i] = 1 / orthonormal(family, q, x[i])[3]
        end
        family ≡ :gh ? (x, w) : ((x .+ 1) ./ 2, w)               # [-1,1] → [0,1]
    end
    reshape(T.(x), :, 1), T.(w)
end

# φ_q(x), φ_q'(x) and Σ_{k<q} φ_k(x)² for the polynomials orthonormal with respect to N(0,1)
# (family :gh) or to the uniform density on [-1,1] (family :le).
function orthonormal(family::Symbol, q::Integer, x::BigFloat)
    a = zero(x); b = one(x); s = zero(x)                         # φ_{k-1}, φ_k
    for k in 0:q-1
        s += b^2
        a, b = b, (family ≡ :gh ? (x * b - sqrt(BigFloat(k)) * a) / sqrt(BigFloat(k + 1)) :
                                  (sqrt(BigFloat((2k + 1) * (2k + 3))) * x * b - k * sqrt(BigFloat(2k + 3) / max(2k - 1, 1)) * a) / (k + 1))
    end
    dφ = family ≡ :gh ? sqrt(BigFloat(q)) * a :
                        q * (x * b - sqrt(BigFloat(2q + 1) / (2q - 1)) * a) / (x^2 - 1)
    b, dφ, s
end

"""
    Quadriceps.extended(family)

The cells of `family` (`:gh` or `:le`) that are stored beyond double precision, and therefore
available as `ghpos(T, d, q)` or `lepos(T, d, q)` for a number type `T` wider than `Float64`: a
vector of named tuples `(family, d, p, n, relerr128)`, sorted by dimension and degree.
`relerr128` is the largest relative monomial error of the rule rounded to IEEE binary128
(`Float128`), measured in much wider arithmetic; its unit is the machine epsilon of that format,
`2^-112 ≈ 1.93e-34`, and a few epsilons is the floor of the format, as a few `1e-16` is for the
`relerr` of the `Float64` rules.
"""
function extended(family::Symbol)
    checkfamily(family)
    sort!([(family = family, d = k[2], p = k[3], n = BIN128[k].n, relerr128 = e) for (k, e) in INDEX128 if k[1] ≡ family && haskey(BIN128, k)];
          by = r -> (r.d, r.p))
end
