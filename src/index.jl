# The catalog of stored rules (data/index.tsv) and the loader for the rule files.

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
* `file` — file name under `data/<family>/`.
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
    file::String
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
                        parse(Float64, r[6]), parse(Float64, r[7]), r[8] == "yes", String(r[9]), String(r[10]))
        index[(info.family, info.d, info.p)] = info
    end
    index
end

function __init__()
    empty!(INDEX)
    merge!(INDEX, readindex(joinpath(datadir(), "index.tsv")))
    nothing
end

# Parse a rule file: one node per line, x_1,…,x_d,w; lines starting with '#' are comments.
function readrule(path::AbstractString)
    rows = Vector{Vector{Float64}}()
    for l in eachline(path)
        s = strip(l)
        (isempty(s) || startswith(s, '#')) && continue
        push!(rows, parse.(Float64, split(s, ',')))
    end
    n = length(rows); d = length(rows[1]) - 1
    X = Matrix{Float64}(undef, n, d); w = Vector{Float64}(undef, n)
    for (i, r) in enumerate(rows)
        length(r) == d + 1 || error("$path: line $i has $(length(r)) fields, expected $(d + 1)")
        X[i, :] .= @view r[1:d]
        w[i] = r[d+1]
    end
    X, w
end

# The stored rule behind a catalog entry, in the normalized frame. The cached arrays are
# shared; callers hand out copies.
function stored(info::RuleInfo)
    key = (info.family, info.d, info.p)
    lock(CACHE_LOCK) do
        get!(CACHE, key) do
            X, w = readrule(joinpath(datadir(), String(info.family), info.file))
            size(X) == (info.n, info.d) || error("$(info.file): size $(size(X)) does not match the catalog")
            X, w
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
