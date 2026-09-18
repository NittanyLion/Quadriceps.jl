# The catalog of stored rules (data/index.tsv) and the reader and writer of data/rules.bin.

"""
    Quadriceps.RuleInfo

Catalog entry of one stored rule. Fields:

* `family` — `:gh` or `:le`;
* `d`, `p`, `n` — dimension, degree of exactness, number of nodes;
* `moller` — Möller's lower bound on `n` for this `(d, p)`, or `-1` where it is not tabulated;
  a rule with `n == moller` is proven minimal;
* `relerr` — largest relative monomial error of the stored `Float64` rule over all monomials
  of total degree `≤ p`, measured when the package data were built
  (see [`Quadriceps.exactness_error`](@ref));
* `minweight` — smallest weight (always `> 0`);
* `interior` — `true` when every node lies inside the integration domain (always `true` for GH);
* `origin` — who the rule belongs to: `"own"`, or `"derived: …"`, `"same-rule: …"`,
  `"transcribed: …"` followed by the published source;
* `source_id` — the same information as the small integer stored in `data/rules.bin`
  (0: own; 3: derived from Diallo and Worku; 10 and up: a published rule; see the format page).
"""
struct RuleInfo
    family::Symbol
    d::Int
    p::Int
    n::Int
    moller::Int
    relerr::Float64
    minweight::Float64
    interior::Bool
    origin::String
    source_id::Int
end

const Catalog = Dict{Tuple{Symbol,Int,Int},RuleInfo}

const FAMILIES = (:gh, :le)
const INDEX = Catalog()
const CACHE = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{Float64},Vector{Float64}}}()
const CACHE_LOCK = ReentrantLock()

datadir() = joinpath(dirname(@__DIR__), "data")

function checkfamily(family::Symbol)
    family ∈ FAMILIES || throw(ArgumentError("unknown family :$family; use :gh or :le"))
    family
end

function readindex(path::AbstractString)
    index = Catalog()
    isfile(path) || return index
    for l in eachline(path)
        (isempty(l) || startswith(l, '#')) && continue
        r = split(l, '\t')
        startswith(r[1], "family") && continue
        info = RuleInfo(Symbol(r[1]), parse(Int, r[2]), parse(Int, r[3]), parse(Int, r[4]), parse(Int, r[5]),
                        parse(Float64, r[6]), parse(Float64, r[7]), r[8] == "yes", String(r[9]), parse(Int, r[10]))
        index[(info.family, info.d, info.p)] = info
    end
    index
end

# --- data/rules.bin: every rule in one binary file ------------------------------------------
# Format QUADRICEPS1 (docs/src/format.md), the GHBEST1 layout of the designed-quadrature
# project with a family field added:
#   header   one ASCII line ending in \n: magic, then key=value tokens
#   index    `cells` records of `index_fields` = 8 little-endian Int64:
#            family (0 = gh, 1 = le), d, p, q, n, offset, nbytes, source_id
#   data     per cell n*(d+1) little-endian Float64, row-major: x_1 … x_d, w for each node
const MAGIC = "QUADRICEPS1"
const NIDX = 8

struct BinEntry
    n::Int
    offset::Int
    nbytes::Int
    source_id::Int
end

const BinIndex = Dict{Tuple{Symbol,Int,Int},BinEntry}
const BIN = BinIndex()

binpath() = joinpath(datadir(), "rules.bin")

function readbinindex(path::AbstractString)
    bin = BinIndex()
    isfile(path) || return bin
    open(path) do io
        header = readuntil(io, '\n')
        tok = split(header)
        (!isempty(tok) && tok[1] == MAGIC) || error("$path: not a $MAGIC file")
        kv = Dict(String.(split(t, '='; limit = 2)) for t in tok[2:end] if occursin('=', t))
        (kv["fmt"] == "1" && kv["endian"] == "little" && kv["float"] == "binary64" && kv["index_fields"] == string(NIDX)) ||
            error("$path: unsupported $MAGIC variant: $header")
        for _ in 1:parse(Int, kv["cells"])
            f, d, p, q, n, off, nb, sid = (Int(ltoh(read(io, Int64))) for _ in 1:NIDX)
            nb == n * (d + 1) * 8 || error("$path: cell d=$d p=$p has nbytes=$nb, expected $(n * (d + 1) * 8)")
            bin[(FAMILIES[f+1], d, p)] = BinEntry(n, off, nb, sid)
        end
    end
    bin
end

# One rule out of the binary file, in the normalized frame.
function readblock(path::AbstractString, e::BinEntry, d::Integer)
    v = Vector{Float64}(undef, e.n * (d + 1))
    open(path) do io
        seek(io, e.offset)
        read!(io, v)
    end
    v .= ltoh.(v)
    A = permutedims(reshape(v, d + 1, e.n))              # row-major on disk
    A[:, 1:d], A[:, d+1]
end

# Write rules (key => (X, w, source_id)) as a QUADRICEPS1 file; cells in (family, d, p) order.
function writebin(path::AbstractString, rules::AbstractDict)
    order = sort!(collect(keys(rules)); by = k -> (findfirst(==(k[1]), FAMILIES), k[2], k[3]))
    header = "$MAGIC fmt=1 endian=little cells=$(length(order)) index_fields=$NIDX float=binary64 " *
             "order=row-major layout=x1..xd,w families=0:gh,1:le\n"
    offset = sizeof(header) + length(order) * NIDX * 8
    open(path, "w") do io
        write(io, header)
        for k in order
            X, _, sid = rules[k]; n, d = size(X); nb = n * (d + 1) * 8
            for v in (findfirst(==(k[1]), FAMILIES) - 1, d, k[3], (k[3] + 1) ÷ 2, n, offset, nb, sid)
                write(io, htol(Int64(v)))
            end
            offset += nb
        end
        for k in order
            X, w, _ = rules[k]
            for i in axes(X, 1)
                for c in axes(X, 2); write(io, htol(Float64(X[i, c]))); end
                write(io, htol(Float64(w[i])))
            end
        end
    end
    path
end

function __init__()
    empty!(INDEX); empty!(BIN); empty!(CACHE)
    merge!(INDEX, readindex(joinpath(datadir(), "index.tsv")))
    merge!(BIN, readbinindex(binpath()))
    for (k, info) in INDEX
        (haskey(BIN, k) && BIN[k].n == info.n) || error("Quadriceps: data/index.tsv and data/rules.bin disagree at $k")
    end
    nothing
end

# The stored rule behind a catalog entry, in the normalized frame. The cached arrays are
# shared; callers hand out copies.
function stored(info::RuleInfo)
    key = (info.family, info.d, info.p)
    lock(CACHE_LOCK) do
        get!(CACHE, key) do
            readblock(binpath(), BIN[key], info.d)
        end
    end
end

"""
    Quadriceps.available(family)

All stored rules of `family` (`:gh` or `:le`), as a vector of [`Quadriceps.RuleInfo`](@ref)
sorted by dimension and degree. Tensor products, which `pragmatic = true` builds on demand,
are not listed.

```julia
[(r.d, r.p, r.n) for r in Quadriceps.available(:gh) if r.d == 3]
```
"""
function available(family::Symbol)
    checkfamily(family)
    sort!([r for r in values(INDEX) if r.family ≡ family]; by = r -> (r.d, r.p))
end
