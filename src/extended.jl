# Rules beyond double precision: the number-type argument of ghpos / lepos, served from two files
# that ship with the package, like data/rules.bin. Nothing is fetched from anywhere.
#
# · data/rules128.bin is a QUADRICEPS1 file with float=binary128: every number is the IEEE
#   binary128 (quadruple precision, 113-bit significand) rounding of the deposit's 80-digit rule,
#   16 little-endian bytes. It is decoded here into 113-bit BigFloats, so the package needs no
#   quadruple-precision dependency; `T.(…)` then gives the caller's type, and for `Float128` of
#   Quadmath.jl that conversion is exact. A type of at most 113 bits (Float128, Double64,
#   BigFloat under `setprecision(BigFloat, 113)`) is served from it.
#
# · data/rules80.bin holds the 80-digit rules themselves (the files rules_extended/<name>.mp.csv of
#   the Zenodo deposit, 10.5281/zenodo.22881864), with float=binary320: the binary128 layout
#   widened to 40 little-endian bytes — 1 sign bit, 15 exponent bits (bias 16383), 304 fraction
#   bits, so a 305-bit significand, about 91 digits, which holds an 80-digit decimal to within
#   2^-305. It is decoded into PREC80-bit BigFloats and converted to the caller's type, so a type
#   wider than 113 bits (BigFloat at its default 256 bits, MultiFloats' Float64x4, …) receives the
#   deposit's 80 digits.
#
# data/index128.tsv is the catalog of both files (they hold the same cells): the deposit file each
# cell was taken from, its SHA-256, and the measured errors relerr128 and relerr80.

const BIN128 = BinIndex()
const BIN80 = BinIndex()
const Ext = @NamedTuple{relerr128::Float64, file::String, sha256::String, relerr80::Float64}
const INDEX128 = Dict{Tuple{Symbol,Int,Int},Ext}()               # key => catalog row of index128.tsv
const CACHE128 = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{BigFloat},Vector{BigFloat}}}()
const CACHE80 = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{BigFloat},Vector{BigFloat}}}()
const BIAS128 = 16383
const PREC80 = 305              # significand bits of binary320, the precision rules80.bin is decoded to
const WORDS320 = 5              # 40 bytes = five UInt64 words, least significant first

bin128path() = joinpath(datadir(), "rules128.bin")
bin80path() = joinpath(datadir(), "rules80.bin")

# --- binary128 ------------------------------------------------------------------------------------

# MPFR operations such as ldexp return a BigFloat at the CURRENT DEFAULT precision, whatever the
# operand's, so every codec below sets the precision it needs explicitly: wide enough to hold the
# operand exactly when encoding, the format's own when decoding.

# The binary128 bit pattern of a BigFloat that holds at most 113 significant bits.
function encode128(x::BigFloat)
    iszero(x) && return UInt128(0)
    e = exponent(x)
    -16382 ≤ e ≤ 16383 || throw(ArgumentError("$x is outside the normal range of binary128"))
    m = setprecision(() -> ldexp(abs(x), 112 - e), BigFloat, max(precision(x), 113))   # in [2^112, 2^113)
    isinteger(m) || throw(ArgumentError("$x has more than 113 significant bits"))
    (UInt128(signbit(x)) << 127) | (UInt128(e + BIAS128) << 112) | (UInt128(BigInt(m)) - (UInt128(1) << 112))
end

function decode128(b::UInt128)
    iszero(b << 1) && return BigFloat(0; precision = 113)
    e = Int((b >> 112) & 0x7fff)
    0 < e < 0x7fff || error("rules128.bin holds a subnormal or non-finite number")
    m = (b & ((UInt128(1) << 112) - 1)) | (UInt128(1) << 112)
    x = setprecision(() -> ldexp(BigFloat(BigInt(m); precision = 113), e - BIAS128 - 112), BigFloat, 113)
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

# --- binary320: the binary128 layout with a 304-bit fraction ---------------------------------------

# The five little-endian words of a BigFloat that holds at most PREC80 significant bits.
function encode320(x::BigFloat)
    b = if iszero(x)
        BigInt(0)
    else
        e = exponent(x)
        -16382 ≤ e ≤ 16383 || throw(ArgumentError("$x is outside the normal range of binary320"))
        m = setprecision(() -> ldexp(abs(x), PREC80 - 1 - e), BigFloat, max(precision(x), PREC80))   # in [2^304, 2^305)
        isinteger(m) || throw(ArgumentError("$x has more than $PREC80 significant bits"))
        (BigInt(signbit(x)) << 319) | (BigInt(e + BIAS128) << 304) | (BigInt(m) - (BigInt(1) << 304))
    end
    ntuple(k -> UInt64((b >> (64 * (k - 1))) & typemax(UInt64)), WORDS320)
end

function decode320(words)
    b = sum(BigInt(words[k]) << (64 * (k - 1)) for k in 1:WORDS320)
    iszero(b & ((BigInt(1) << 319) - 1)) && return BigFloat(0; precision = PREC80)
    e = Int((b >> 304) & 0x7fff)
    0 < e < 0x7fff || error("rules80.bin holds a subnormal or non-finite number")
    m = (b & ((BigInt(1) << 304) - 1)) | (BigInt(1) << 304)
    x = setprecision(() -> ldexp(BigFloat(m; precision = PREC80), e - BIAS128 - 304), BigFloat, PREC80)
    iszero(b >> 319) ? x : -x
end

function readblock80(path::AbstractString, e::BinEntry, d::Integer)
    v = Vector{UInt64}(undef, e.n * (d + 1) * WORDS320)
    open(path) do io
        seek(io, e.offset)
        read!(io, v)
    end
    v .= ltoh.(v)
    vals = [decode320(view(v, WORDS320 * (j - 1) + 1:WORDS320 * j)) for j in 1:e.n*(d+1)]
    A = permutedims(reshape(vals, d + 1, e.n))                   # row-major on disk
    A[:, 1:d], A[:, d+1]
end

# --- the catalog and the two stores -----------------------------------------------------------------

# index128.tsv: family d p n relerr128 extendedfile sha256 relerr80 (the last three were added before the first registered version).
function readindex128(path::AbstractString)
    out = Dict{Tuple{Symbol,Int,Int},Ext}()
    isfile(path) || return out
    for l in eachline(path)
        (isempty(l) || startswith(l, '#') || startswith(l, "family")) && continue
        r = split(l, '\t')
        length(r) ≥ 8 || error("Quadriceps: data/index128.tsv has a row with $(length(r)) fields; 8 are needed")
        out[(Symbol(r[1]), parse(Int, r[2]), parse(Int, r[3]))] =
            (relerr128 = parse(Float64, r[5]), file = String(r[6]), sha256 = String(r[7]), relerr80 = parse(Float64, r[8]))
    end
    out
end

# Since the 2026-09-19 QUAD_ONLY policy the catalog holds no cell without a binary128 rule (the
# builder leaves such a cell out), so the errors below guard against damaged data, not a case a
# caller can reach through ghpos / lepos.
function checkwide(bin::BinIndex, info::RuleInfo, what::AbstractString)
    key = (info.family, info.d, info.p)
    (haskey(bin, key) && bin[key].n == info.n) ||
        throw(ArgumentError("the rule for family :$(info.family), d = $(info.d), p = $(info.p) has no $what; " *
                            "see Quadriceps.extended(:$(info.family)) for the cells held beyond double precision"))
    key
end

# The stored rule in quadruple precision (113-bit BigFloats, shared: callers convert, which copies).
function stored128(info::RuleInfo)
    key = checkwide(BIN128, info, "binary128 rule")
    lock(CACHE_LOCK) do
        get!(CACHE128, key) do
            readblock128(bin128path(), BIN128[key], info.d)
        end
    end
end

# The stored rule to 80 digits (PREC80-bit BigFloats, shared: callers convert, which copies).
function stored80(info::RuleInfo)
    key = checkwide(BIN80, info, "80-digit rule")
    lock(CACHE_LOCK) do
        get!(CACHE80, key) do
            readblock80(bin80path(), BIN80[key], info.d)
        end
    end
end

# A number type that Float64 data serve in full.
narrow(::Type{T}) where {T} = T ≡ Float64 || T ≡ Float32 || T ≡ Float16

# Significant bits of T; a type without a `precision` method is taken to be wide. For BigFloat this
# is the current default precision, so `setprecision(BigFloat, 113)` selects the binary128 data.
bits(::Type{T}) where {T} = hasmethod(precision, Tuple{Type{T}}) ? precision(T) : typemax(Int)
wide(::Type{T}) where {T} = bits(T) > 113

function stored(::Type{T}, info::RuleInfo) where {T<:AbstractFloat}
    X, w = narrow(T) ? stored(info) : wide(T) ? stored80(info) : stored128(info)
    T.(X), T.(w)
end

# One-dimensional Gauss rules in T: the Float64 nodes of FastGaussQuadrature.jl, refined by
# Newton's method on the orthonormal polynomial in 256-bit arithmetic (PREC80 bits for a wide T);
# weights from the Christoffel function, 1 / Σ_{k<q} φ_k(x)², which sum to 1 for a probability weight.
function gauss1d(::Type{T}, family::Symbol, q::Integer) where {T<:AbstractFloat}
    narrow(T) && return (r = gauss1d(family, q); (T.(r[1]), T.(r[2])))
    x, w = setprecision(BigFloat, wide(T) ? PREC80 : 256) do
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

The cells of `family` (`:gh` or `:le`) that are stored beyond double precision — every stored
cell — and so available as `ghpos(T, d, q)` or `lepos(T, d, q)` for a number type `T` wider than
`Float64`: a vector of named tuples `(family, d, p, n, relerr128, relerr80)`, sorted by dimension
and degree. Both errors are the largest relative monomial error of the rule, measured in much
wider arithmetic: `relerr128` of the rule rounded to IEEE binary128 (`Float128`; the unit is the
machine epsilon of that format, `2^-112 ≈ 1.93e-34`, and a few epsilons is the floor of the
format, as a few `1e-16` is for the `relerr` of the `Float64` rules), `relerr80` of the 80-digit
rule that a type wider than 113 bits receives, as the deposit measured it.
"""
function extended(family::Symbol)
    checkfamily(family)
    sort!([(family = family, d = k[2], p = k[3], n = BIN128[k].n, relerr128 = e.relerr128, relerr80 = e.relerr80)
           for (k, e) in INDEX128 if k[1] ≡ family && haskey(BIN128, k)];
          by = r -> (r.d, r.p))
end
