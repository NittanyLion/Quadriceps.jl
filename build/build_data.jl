# build_data.jl — fill data/ from the rule bank of the designed-quadrature project.
#
#   julia --project=build build/build_data.jl [--force] [path to the project's sync folder]
#
# For every (family, d, p) the bank holds (family = hermite → gh, legendre → le), take the
# smallest rule that passes the checks below, and keep it if it is worth storing:
#
#   checks   all weights > 0; largest relative monomial error (Quadriceps.exactness_error,
#            Float64) below GATE; for Le, all nodes inside [0,1]^d is recorded, not required.
#   keep     n is below the cheapest tensor product of lower-dimensional kept rules and
#            one-dimensional Gauss rules, or n equals Möller's bound (proven minimal).
#            The bank's own tensor fills therefore stay out: `pragmatic = true` rebuilds them.
#
# Origin of each rule: rules/literature.tsv (transcribed / same-rule) and the parent chains in
# rules/lineage/*.tsv (a chain that ends at a dw: file is derived from Diallo and Worku's rule).
#
# The new data are built in data.new/ and replace data/ only if they are no worse than what
# data/index.tsv holds now. Exit codes (build/update.sh, which runs this unattended, relies on
# them; --force overrides 3 and 4):
#
#   0   data/ is current: unchanged, or replaced (one "CHANGE: …" line per cell that changed)
#   3   regression — a stored cell would disappear or get more nodes (a bank file caught
#       half-synced, or a rule withdrawn from the bank); data/ is left alone
#   4   the set of rules that are not the author's own would change; NOTICE.md and
#       docs/src/credits.md need a person, so data/ is left alone
#
# The rules go into ONE binary file, data/rules.bin (format QUADRICEPS1, docs/src/format.md);
# data/index.tsv is the catalog that goes with it. No per-rule text files are written.
#
# Also rewrites docs/src/rules.md and the generated blocks of README.md. The output carries no
# dates, so a rebuild from an unchanged bank changes no file.

using Quadriceps, Printf, SHA
using Quadriceps: RuleInfo, Catalog, cheapest, readindex, readbinindex, readblock, writebin, exactness_error

const GATE = 1e-11
const PKG  = dirname(@__DIR__)
const FORCE = "--force" ∈ ARGS
const POS  = filter(≠("--force"), ARGS)
const SYNC = length(POS) ≥ 1 ? POS[1] : joinpath(homedir(), "Dropbox", "oldDesignedQuadrature-sync")
const NEW  = joinpath(PKG, "data.new")
const BANK = joinpath(SYNC, "project", "julia", "rules")
const FAM  = Dict("hermite" => :gh, "legendre" => :le)

# A bank file: one node per line, x_1,…,x_d,w; lines starting with '#' are comments.
function readrule(path::AbstractString)
    M = [parse.(Float64, split(strip(l), ',')) for l in eachline(path) if !isempty(strip(l)) && !startswith(strip(l), '#')]
    d = length(M[1]) - 1
    all(r -> length(r) == d + 1, M) || error("$path: ragged rows")
    A = permutedims(reduce(hcat, M))
    A[:, 1:d], A[:, d+1]
end

# source_id, the integer form of the origin stored in rules.bin. 0–15 as in GHBEST1
# (bestknown_gh/NOTES.md §7); 3 and 16 are new here. Keep docs/src/format.md in step.
const SOURCE_IDS = ["Stroud and Secrest 1963" => 11, "Stroud 1971" => 10, "Haegemans and Piessens 1976" => 12,
                    "Haegemans and Piessens 1977" => 13, "Cools and Haegemans 1988" => 14, "Konyaev 1977" => 15,
                    "Festa and Sommariva 2012" => 16]
function source_id(org)
    org == "own" && return 0
    startswith(org, "derived: Diallo and Worku") && return 3
    for (name, id) in SOURCE_IDS
        occursin(name, org) && return id
    end
    99                                                    # a published source without an id yet
end

rows(path) = [split(l, '\t') for l in eachline(path) if !isempty(l) && !startswith(l, '#')]

# --- bank: every candidate per cell, smallest first ---------------------------------------
cands = Dict{Tuple{Symbol,Int,Int},Vector{Tuple{Int,String}}}()
for f in readdir(BANK)
    m = match(r"^(hermite|legendre)_d(\d+)_p(\d+)_n(\d+)\.csv$", f)
    m ≡ nothing && continue
    push!(get!(cands, (FAM[m[1]], parse(Int, m[2]), parse(Int, m[3])), Tuple{Int,String}[]), (parse(Int, m[4]), f))
end
foreach(sort!, values(cands))

moller = Dict((parse(Int, r[1]), parse(Int, r[2])) => parse(Int, r[3]) for r in rows(joinpath(BANK, "..", "symq", "moller.tsv")))

# --- origin ---------------------------------------------------------------------------------
register = Dict(String(r[1]) => (String(r[2]), String(r[3])) for r in rows(joinpath(BANK, "literature.tsv")) if length(r) ≥ 3)
parent = Dict{String,String}()
for f in readdir(joinpath(BANK, "lineage"); join = true)
    endswith(f, ".tsv") || continue
    for r in rows(f)
        length(r) ≥ 2 && r[2] ≠ "none" && r[1] ≠ r[2] && (parent[String(r[1])] = String(r[2]))
    end
end
function dwroot(f)
    seen = Set{String}()
    while haskey(parent, f) && f ∉ seen
        push!(seen, f); f = parent[f]
        startswith(f, "dw:") && return replace(f[4:end], ".csv" => "")
    end
    nothing
end
function origin(f)
    haskey(register, f) && return "$(register[f][1]): $(register[f][2])"
    dw = dwroot(f)
    dw ≡ nothing ? "own" : "derived: Diallo and Worku 2026, elimination started from their rule $dw for the same cell"
end

# --- select ---------------------------------------------------------------------------------
const WEIGHT = Dict(:gh => "standard normal N(0, I_d), density (2π)^(-d/2) exp(-|x|²/2) on R^d; weights sum to 1",
                    :le => "uniform density on [0,1]^d; weights sum to 1")
index = Catalog()
dropped = String[]
rm(NEW; recursive = true, force = true); mkpath(NEW)
rules = Dict{Tuple{Symbol,Int,Int},Tuple{Matrix{Float64},Vector{Float64},Int}}()
provenance = Dict{Tuple{Symbol,Int,Int},Tuple{String,String}}()
for key in sort!(collect(keys(cands)); by = k -> (k[2], k[3], k[1]))      # lower dimensions first
    fam, d, p = key
    for (n, f) in cands[key]
        X, w = readrule(joinpath(BANK, f))
        size(X) == (n, d) || (push!(dropped, "$f: size $(size(X)) does not match its name"); continue)
        minimum(w) > 0 || (push!(dropped, "$f: nonpositive weight $(minimum(w))"); continue)
        err = exactness_error(X, w, p, fam)
        err < GATE || (push!(dropped, @sprintf("%s: relative error %.1e", f, err)); continue)
        mb = get(moller, (d, p), -1)
        tcost = d == 1 ? typemax(Int) : minimum(cheapest(index, fam, j, p)[2] * cheapest(index, fam, d - j, p)[2] for j in 1:d÷2)
        if !(n < tcost || n == mb)
            push!(dropped, "$f: not below the tensor product of lower-dimensional rules ($tcost nodes)")
            break
        end
        inside = fam ≡ :gh || all(x -> 0 < x < 1, X)
        org = origin(f)
        rules[key] = (X, w, source_id(org))
        provenance[key] = (f, bytes2hex(sha256(read(joinpath(BANK, f)))))
        index[key] = RuleInfo(fam, d, p, n, mb, err, minimum(w), inside, org, source_id(org))
        break
    end
end

# --- index ----------------------------------------------------------------------------------
infos = sort!(collect(values(index)); by = r -> (r.family, r.d, r.p))
open(joinpath(NEW, "index.tsv"), "w") do io
    println(io, "# Quadriceps.jl rule catalog — written by build/build_data.jl; do not edit by hand")
    println(io, "# the rules themselves are in rules.bin; bankfile and sha256 name the bank file each rule was taken from")
    println(io, "family\td\tp\tn\tmoller\trelerr\tminweight\tinterior\torigin\tsource_id\tbankfile\tsha256")
    for r in infos
        @printf(io, "%s\t%d\t%d\t%d\t%d\t%.3e\t%.6e\t%s\t%s\t%d\t%s\t%s\n", r.family, r.d, r.p, r.n, r.moller, r.relerr,
                r.minweight, r.interior ? "yes" : "no", r.origin, r.source_id, provenance[(r.family, r.d, r.p)]...)
    end
end
writebin(joinpath(NEW, "rules.bin"), rules)

# --- compare with what is stored now, then swap -----------------------------------------------
old = readindex(joinpath(PKG, "data", "index.tsv"))
kind(r) = first(split(r.origin, ':'))
cellname(k) = "$(k[1] ≡ :gh ? "GH" : "Le") d=$(k[2]) p=$(k[3])"
oldbin = readbinindex(joinpath(PKG, "data", "rules.bin"))
same(k) = haskey(oldbin, k) && readblock(joinpath(PKG, "data", "rules.bin"), oldbin[k], k[2]) == rules[k][1:2]
worse = String[]; credit = String[]; changes = String[]
for (k, o) in old
    haskey(index, k) || (push!(worse, "$(cellname(k)): n=$(o.n) would disappear"); continue)
    index[k].n > o.n && push!(worse, "$(cellname(k)): n=$(o.n) would become $(index[k].n)")
end
for k in sort!(collect(union(keys(old), keys(index))); by = k -> (k[1], k[2], k[3]))
    ko = haskey(old, k) ? kind(old[k]) : "own"; kn = haskey(index, k) ? kind(index[k]) : "own"
    ko ≠ kn && push!(credit, "$(cellname(k)): $ko → $kn")
    haskey(index, k) || continue
    r = index[k]
    if !haskey(old, k)
        push!(changes, "$(cellname(k)): new, n=$(r.n)")
    elseif old[k].n ≠ r.n
        push!(changes, "$(cellname(k)): n $(old[k].n) → $(r.n)")
    elseif !same(k)
        push!(changes, "$(cellname(k)): n=$(r.n), rule replaced")
    end
end
fills = count(contains("not below the tensor product"), dropped)
fills > 0 && println("left out: $fills bank rules that a tensor product of lower-dimensional rules matches")
filter!(!contains("not below the tensor product"), dropped)
isempty(dropped) || println("left out:\n  ", join(dropped, "\n  "))
if !FORCE && !isempty(old) && !isempty(worse)
    println("REGRESSION, data/ left alone:\n  ", join(worse, "\n  ")); rm(NEW; recursive = true); exit(3)
end
if !FORCE && !isempty(old) && !isempty(credit)
    println("CREDIT CHANGE, data/ left alone (update NOTICE.md and docs/src/credits.md, then --force):\n  ", join(credit, "\n  "))
    rm(NEW; recursive = true); exit(4)
end
rm(joinpath(PKG, "data"); recursive = true, force = true)
mv(NEW, joinpath(PKG, "data"))
foreach(c -> println("CHANGE: ", c), changes)

# --- docs/src/rules.md ------------------------------------------------------------------------
open(joinpath(PKG, "docs", "src", "rules.md"), "w") do io
    println(io, "# Stored rules\n")
    println(io, "Written by `build/build_data.jl`. `q`: the argument of `ghpos(d, q)` and `lepos(d, q)`;")
    println(io, "`p = 2q - 1`: degree of exactness; `n`: number of nodes; `ρ = n^(1/d) / q`: the node count")
    println(io, "relative to the `q^d` Gauss product grid (1.00 is that grid, smaller")
    println(io, "is better); Möller: Möller's lower bound on `n` (**bold** `n`: bound attained, proven minimal);")
    println(io, "rel. err.: largest relative monomial error of the stored `Float64` rule; origin: `own`, or the")
    println(io, "published rule it is, or descends from (see [Credits](credits.md)).\n")
    for (fam, title) in ((:gh, "GH — Gaussian weight"), (:le, "Le — uniform weight on the cube"))
        println(io, "## $title\n")
        for d in sort!(unique(r.d for r in infos if r.family ≡ fam))
            println(io, "### d = $d\n\n| q | p | n | ρ | Möller | rel. err. | origin |\n|---:|---:|---:|---:|---:|---:|:---|")
            for r in infos
                (r.family ≡ fam && r.d == d) || continue
                ρ = r.n^(1 / d) / ((r.p + 1) ÷ 2)
                @printf(io, "| %d | %d | %s | %.2f | %s | %.1e | %s |\n", (r.p + 1) ÷ 2, r.p, r.n == r.moller ? "**$(r.n)**" : string(r.n), ρ,
                        r.moller < 0 ? "—" : string(r.moller), r.relerr, r.origin)
            end
            println(io)
        end
    end
end

# --- generated blocks of README.md ------------------------------------------------------------
# <!-- BEGIN GENERATED name --> … <!-- END GENERATED name -->; everything between is replaced.
function regenerate(text, name, body)
    b = "<!-- BEGIN GENERATED $name -->"; e = "<!-- END GENERATED $name -->"
    i = findfirst(b, text); j = findfirst(e, text)
    (i ≡ nothing || j ≡ nothing) && error("README.md has no generated block '$name'")
    text[1:last(i)] * "\n" * body * "\n" * text[first(j):end]
end
top(fam, d) = (rs = [r for r in infos if r.family ≡ fam && r.d == d]; (maximum(r.p for r in rs), maximum(r.n for r in rs)))
nn(fam, d, p) = cheapest(index, fam, d, p)[2]
coverage = let io = IOBuffer()
    println(io, "$(length(infos)) rules are stored: for example $(nn(:gh, 5, 9)) nodes instead of $(5^5) for the Gaussian weight at")
    println(io, "`d = 5, q = 5`, and $(nn(:le, 5, 21)) instead of $(11^5) for the cube at `d = 5, q = 11`.\n")
    println(io, "| `d` | GH | largest GH rule | Le | largest Le rule |\n|---|---|---|---|---|")
    for d in sort!(unique(r.d for r in infos))
        (pg, ng), (pl, nl) = top(:gh, d), top(:le, d)
        println(io, "| $d | `q ≤ $((pg + 1) ÷ 2)` (`p ≤ $pg`) | $ng nodes | `q ≤ $((pl + 1) ÷ 2)` (`p ≤ $pl`) | $nl nodes |")
    end
    String(take!(io))
end
pragmatic = """
```julia
ghpos(7, 5)                             # ArgumentError: no stored rule in seven dimensions
X, w = ghpos(7, 5; pragmatic = true)    # $(nn(:gh, 7, 9)) nodes; the product grid has $(5^7)
X, w = lepos(3, 24; pragmatic = true)   # $(nn(:le, 3, 47)) nodes; the product grid has $(24^3)

Quadriceps.nnodes(:gh, 10, 3; pragmatic = true)     # $(nn(:gh, 10, 5)), without building the rule
Quadriceps.ruleinfo(:gh, 7, 5; pragmatic = true)    # the factors, with their origins
```"""
count(k) = sum(startswith(r.origin, k) for r in infos)
credits = "$(count("own")) of the $(length(infos)) rules were computed by the author. $(count("transcribed") + count("same-rule")) are rules from the literature (copied in,\n" *
          "or found again by the author's search and recognized), and $(count("derived")) Le rules were obtained by node\n" *
          "elimination started from Diallo and Worku's published rules."
readme = read(joinpath(PKG, "README.md"), String)
updated = regenerate(regenerate(regenerate(readme, "coverage", coverage), "pragmatic", pragmatic), "credits", credits)
updated ≠ readme && write(joinpath(PKG, "README.md"), updated)

println("$(length(infos)) rules stored; ", sum(r.n for r in infos), " nodes; ", isempty(changes) ? "no change" : "$(length(changes)) cells changed")
