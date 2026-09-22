# Rules beyond double precision: the number-type argument of ghpos / lepos, served from two sources.
#
# · data/rules128.bin is a QUADRICEPS1 file with float=binary128: every number is the IEEE
#   binary128 (quadruple precision, 113-bit significand) rounding of the deposit's 80-digit rule,
#   16 little-endian bytes. It is decoded here into 113-bit BigFloats, so the package needs no
#   quadruple-precision dependency; `T.(…)` then gives the caller's type, and for `Float128` of
#   Quadmath.jl that conversion is exact. data/index128.tsv is its catalog. A type of at most 113
#   bits (Float128, Double64, BigFloat under `setprecision(BigFloat, 113)`) is served from it.
#
# · The 80-digit rules themselves are the files rules_extended/<name>.mp.csv of the Zenodo deposit
#   (record 10.5281/zenodo.22881864, archives gh.tar.xz and le.tar.xz), declared in Artifacts.toml
#   as lazy artifacts: nothing is downloaded until a caller asks for a type wider than 113 bits
#   (BigFloat at its default 256 bits, MultiFloats' Float64x4, …), then the family's archive is
#   fetched once, verified by Pkg against its SHA-256 and tree hash, and kept in the artifact store.
#   The files are parsed into PREC80-bit BigFloats and converted to the caller's type, so such a
#   result carries the deposit's 80 digits. data/index128.tsv names the file of each cell, its
#   SHA-256 and the deposit's measured error of the 80-digit rule (relerr80).

const BIN128 = BinIndex()
const Ext = @NamedTuple{relerr128::Float64, file::String, sha256::String, relerr80::Float64}
const INDEX128 = Dict{Tuple{Symbol,Int,Int},Ext}()               # key => catalog row of index128.tsv
const CACHE128 = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{BigFloat},Vector{BigFloat}}}()
const CACHE80 = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{BigFloat},Vector{BigFloat}}}()
const BIAS128 = 16383
const PREC80 = 320              # bits held per number of an 80-digit file (80 digits ≈ 266 bits)
const DEPOSIT_DOI = "10.5281/zenodo.22881864"

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

# index128.tsv: family d p n relerr128 extendedfile sha256 relerr80 (the last three since v0.2).
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

# The stored rule in quadruple precision (113-bit BigFloats, shared: callers convert, which copies).
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

# --- 80 digits: the deposit's files, as lazy artifacts -------------------------------------------

# The unpacked archive of a family: gh/ or le/ inside the artifact directory. The first call for a
# family downloads its archive from Zenodo; Pkg verifies it. Failure (no network, a blocked
# download) is reported with the way around it.
function depositdir(family::Symbol)
    root = try
        family ≡ :gh ? artifact"quadriceps_gh" : artifact"quadriceps_le"
    catch e
        error("Quadriceps: the 80-digit $(family ≡ :gh ? "GH" : "Le") rules are fetched once from the Zenodo deposit " *
              "(record $DEPOSIT_DOI, $(family ≡ :gh ? "2.4" : "15") MB) and the download failed: $(sprint(showerror, e)). " *
              "Without network access, ask for a type of at most 113 bits — Float128, or BigFloat under " *
              "setprecision(BigFloat, 113) — which the package serves from its own data.")
    end
    joinpath(root, family ≡ :gh ? "gh" : "le")
end

# One deposit file, rules_extended/<name>.mp.csv: comment lines, a header line x1,…,xd,w, then n
# rows of d + 1 decimal numbers with 80 significant digits, in the row order of the stored rule.
function readextended(path::AbstractString, d::Integer, n::Integer)
    isfile(path) || error("Quadriceps: $path is missing from the deposit archive")
    X = Matrix{BigFloat}(undef, n, d); w = Vector{BigFloat}(undef, n); i = 0
    for l in eachline(path)
        occursin(r"^\s*[-+0-9.]", l) || continue
        s = split(strip(l), ',')
        length(s) == d + 1 || error("Quadriceps: $path has a row with $(length(s)) fields, not $(d + 1)")
        i += 1
        i ≤ n || error("Quadriceps: $path has more than $n rows")
        for k in 1:d
            X[i, k] = BigFloat(String(s[k]); precision = PREC80)
        end
        w[i] = BigFloat(String(s[d+1]); precision = PREC80)
    end
    i == n || error("Quadriceps: $path has $i rows, not $n")
    X, w
end

# The stored rule to 80 digits (PREC80-bit BigFloats, shared: callers convert, which copies).
function stored80(info::RuleInfo)
    key = (info.family, info.d, info.p)
    haskey(INDEX128, key) && BIN128[key].n == info.n ||
        throw(ArgumentError("the rule for family :$(info.family), d = $(info.d), p = $(info.p) has no extended-precision file"))
    lock(CACHE_LOCK) do
        get!(CACHE80, key) do
            readextended(joinpath(depositdir(info.family), "rules_extended", INDEX128[key].file), info.d, info.n)
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
