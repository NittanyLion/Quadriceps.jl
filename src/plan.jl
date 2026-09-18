# Which rule answers a request (d, p): a stored one, or the cheapest tensor product.

# One factor of a plan: a stored rule (info ≢ nothing) or the q-node one-dimensional Gauss rule.
struct Atom
    d::Int
    n::Int
    info::Union{RuleInfo,Nothing}
end

# Nodes of the one-dimensional Gauss rule of degree ≥ p.
gaussq(p::Integer) = p ÷ 2 + 1

# The smallest stored rule of dimension d whose degree is at least p (ties: the lowest degree),
# or nothing. A rule exact to degree p′ ≥ p is a rule of degree p.
function beststored(index::Catalog, family::Symbol, d::Integer, p::Integer)
    best = nothing
    for info in values(index)
        (info.family ≡ family && info.d == d && info.p ≥ p) || continue
        if best ≡ nothing || info.n < best.n || (info.n == best.n && info.p < best.p)
            best = info
        end
    end
    best
end

# The single-rule answer for (d, p): Gauss in one dimension, a stored rule otherwise.
function atom(index::Catalog, family::Symbol, d::Integer, p::Integer)
    d == 1 && return Atom(1, gaussq(p), nothing)
    info = beststored(index, family, d, p)
    info ≡ nothing ? nothing : Atom(d, info.n, info)
end

# The cheapest way to cover d dimensions at degree p with a product of atoms, by dynamic
# programming over the dimension (a product of products is a product, so splitting in two
# is enough). A single stored rule wins ties. Node counts are BigInt: q^d overflows early.
function cheapest(index::Catalog, family::Symbol, d::Integer, p::Integer)
    cost = Vector{BigInt}(undef, d)
    parts = Vector{Vector{Atom}}(undef, d)
    for k in 1:d
        a = atom(index, family, k, p)
        found = a ≢ nothing
        found && (cost[k] = a.n; parts[k] = [a])
        for j in 1:k÷2
            c = cost[j] * cost[k-j]
            if !found || c < cost[k]
                cost[k] = c; parts[k] = vcat(parts[j], parts[k-j]); found = true
            end
        end
    end
    parts[d], cost[d]
end

function norule(family::Symbol, d::Integer, p::Integer, index::Catalog)
    ps = [info.p for info in values(index) if info.family ≡ family && info.d == d]
    have = isempty(ps) ? "no rules are stored for d = $d" :
           "stored rules for d = $d reach q = $((maximum(ps) + 1) ÷ 2), p = $(maximum(ps))"
    at = isodd(p) ? "q = $((p + 1) ÷ 2) (p = $p)" : "p = $p"
    ArgumentError("no stored positive-weight $(family ≡ :gh ? "GH" : "Le") rule for d = $d, $at: $have; " *
                  "pragmatic = true returns the cheapest tensor product of lower-dimensional rules instead")
end

function plan(index::Catalog, family::Symbol, d::Integer, p::Integer, pragmatic::Bool)
    checkfamily(family)
    d ≥ 1 || throw(ArgumentError("the dimension d must be at least 1, got $d"))
    p ≥ 0 || throw(ArgumentError("the degree p must be nonnegative, got $p"))
    if pragmatic
        cheapest(index, family, d, p)
    else
        a = atom(index, family, d, p)
        a ≡ nothing && throw(norule(family, d, p, index))
        [a], BigInt(a.n)
    end
end

# The one-dimensional Gauss rule in the normalized frame: N(0,1) for GH, uniform on [0,1] for Le.
function gauss1d(family::Symbol, q::Integer)
    if family ≡ :gh
        x, w = gausshermite(q; normalize = true)
        reshape(x, :, 1), w
    else
        x, w = gausslegendre(q)
        reshape((x .+ 1) ./ 2, :, 1), w ./ 2
    end
end

# Tensor product of rules; the first factor varies slowest.
function tensor(rules)
    n = prod(length(w) for (_, w) in rules)
    d = sum(size(X, 2) for (X, _) in rules)
    X = Matrix{Float64}(undef, n, d); w = ones(n)
    rep = n; col = 0
    for (Xk, wk) in rules
        nk = length(wk); rep ÷= nk
        for i in 1:n
            w[i] *= wk[((i - 1) ÷ rep) % nk + 1]
        end
        for c in axes(Xk, 2), i in 1:n
            X[i, col+c] = Xk[((i - 1) ÷ rep) % nk + 1, c]
        end
        col += size(Xk, 2)
    end
    X, w
end

function materialize(family::Symbol, parts::Vector{Atom}, n::BigInt)
    d = sum(a.d for a in parts)
    bytes = n * (d + 1) * 8
    bytes ≤ Sys.total_memory() ||
        throw(ArgumentError("the cheapest rule for this request has $n nodes, which does not fit in memory"))
    rules = [a.info ≡ nothing ? gauss1d(family, a.n) : stored(a.info) for a in parts]
    length(rules) == 1 ? (copy(rules[1][1]), copy(rules[1][2])) : tensor(rules)
end
