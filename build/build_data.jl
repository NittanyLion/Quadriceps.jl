# build_data.jl — fill data/ from the rule bank of the designed-quadrature project.
#
#   julia --project=build build/build_data.jl [path to the project's sync folder]
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
# Writes data/<family>/*.csv, data/index.tsv and docs/src/rules.md.

using Quadriceps, Printf, SHA, Dates
using Quadriceps: RuleInfo, Catalog, cheapest, readrule, exactness_error

const GATE = 1e-11
const PKG  = dirname(@__DIR__)
const SYNC = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(homedir(), "Dropbox", "oldDesignedQuadrature-sync")
const BANK = joinpath(SYNC, "project", "julia", "rules")
const FAM  = Dict("hermite" => :gh, "legendre" => :le)

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
for fam in (:gh, :le)
    dir = joinpath(PKG, "data", String(fam)); rm(dir; recursive = true, force = true); mkpath(dir)
end
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
        name = "$(fam)_d$(d)_p$(p)_n$(n).csv"
        org = origin(f)
        open(joinpath(PKG, "data", String(fam), name), "w") do io
            println(io, "# Quadriceps.jl — $(fam ≡ :gh ? "GH" : "Le") rule, d = $d, degree p = $p, n = $n nodes, all weights positive")
            println(io, "# weight: ", WEIGHT[fam])
            println(io, "# columns: ", join(["x$k" for k in 1:d], ","), ",w")
            println(io, "# origin: ", org)
            println(io, "# source: bank file $f, sha256 ", bytes2hex(sha256(read(joinpath(BANK, f)))))
            for i in 1:n
                println(io, join((repr(X[i, k]) for k in 1:d), ","), ",", repr(w[i]))
            end
        end
        index[key] = RuleInfo(fam, d, p, n, mb, err, minimum(w), inside, org, name)
        @printf("%s d=%d p=%-2d n=%-6d err %.1e  %s\n", fam, d, p, n, err, first(org, 40))
        break
    end
end

# --- index ----------------------------------------------------------------------------------
infos = sort!(collect(values(index)); by = r -> (r.family, r.d, r.p))
open(joinpath(PKG, "data", "index.tsv"), "w") do io
    println(io, "# Quadriceps.jl rule catalog — written by build/build_data.jl on $(today()); do not edit by hand")
    println(io, "family\td\tp\tn\tmoller\trelerr\tminweight\tinterior\torigin\tfile")
    for r in infos
        @printf(io, "%s\t%d\t%d\t%d\t%d\t%.3e\t%.6e\t%s\t%s\t%s\n", r.family, r.d, r.p, r.n, r.moller, r.relerr,
                r.minweight, r.interior ? "yes" : "no", r.origin, r.file)
    end
end

# --- docs/src/rules.md ------------------------------------------------------------------------
open(joinpath(PKG, "docs", "src", "rules.md"), "w") do io
    println(io, "# Stored rules\n")
    println(io, "Written by `build/build_data.jl` on $(today()). `n`: number of nodes; `ρ = n^(1/d) / q` with")
    println(io, "`q = (p+1)/2`: the node count relative to the `q^d` Gauss product grid (1.00 is that grid, smaller")
    println(io, "is better); Möller: Möller's lower bound on `n` (**bold** `n`: bound attained, proven minimal);")
    println(io, "rel. err.: largest relative monomial error of the stored `Float64` rule; origin: `own`, or the")
    println(io, "published rule it is, or descends from (see [Credits](credits.md)).\n")
    for (fam, title) in ((:gh, "GH — Gaussian weight"), (:le, "Le — uniform weight on the cube"))
        println(io, "## $title\n")
        for d in sort!(unique(r.d for r in infos if r.family ≡ fam))
            println(io, "### d = $d\n\n| p | n | ρ | Möller | rel. err. | origin |\n|---:|---:|---:|---:|---:|:---|")
            for r in infos
                (r.family ≡ fam && r.d == d) || continue
                ρ = r.n^(1 / d) / ((r.p + 1) ÷ 2)
                @printf(io, "| %d | %s | %.2f | %s | %.1e | %s |\n", r.p, r.n == r.moller ? "**$(r.n)**" : string(r.n), ρ,
                        r.moller < 0 ? "—" : string(r.moller), r.relerr, r.origin)
            end
            println(io)
        end
    end
end

println("\n$(length(infos)) rules stored; ", sum(r.n for r in infos), " nodes")
isempty(dropped) || println("\nleft out:\n  ", join(dropped, "\n  "))
