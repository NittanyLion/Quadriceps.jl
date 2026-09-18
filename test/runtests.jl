using Quadriceps, Test
using FastGaussQuadrature: gausshermite, gausslegendre
using Quadriceps: available, nnodes, ruleinfo, exactness_error

# Every stored rule passed a 1e-11 gate when the data were built; the catalog records the error.
const GATE = 1e-11

@testset "Quadriceps" begin

    @testset "catalog" begin
        @test !isempty(available(:gh)) && !isempty(available(:le))
        @test all(r -> isodd(r.p) && r.d ≥ 2 && r.n ≥ 1, vcat(available(:gh), available(:le)))
        @test_throws ArgumentError available(:la)
    end

    @testset "stored $fam rules" for (fam, f) in ((:gh, ghpos), (:le, lepos))
        for r in available(fam)
            X, w = f(r.d, r.p)
            @test size(X) == (r.n, r.d) && length(w) == r.n
            @test all(>(0), w)
            @test sum(w) ≈ 1
            @test exactness_error(X, w, r.p, fam) < GATE
            fam ≡ :le && @test all(x -> 0 < x < 1, X)
            @test nnodes(fam, r.d, r.p) ≤ r.n               # a higher degree may be cheaper, never dearer
            @test nnodes(fam, r.d, r.p; pragmatic = true) == nnodes(fam, r.d, r.p)   # no product beats a stored rule
        end
    end

    @testset "one dimension" begin
        x, v = gausshermite(4; normalize = true)
        X, w = ghpos(1, 7)
        @test size(X) == (4, 1) && vec(X) == x && w == v
        x, v = gausshermite(4)
        X, w = ghpos(1, 6; normalize = false)               # even degree: next odd
        @test vec(X) ≈ x && w ≈ v
        x, v = gausslegendre(5)
        X, w = lepos(1, 9; normalize = false)
        @test vec(X) ≈ x && w ≈ v
        X, w = lepos(1, 9)
        @test sum(w) ≈ 1 && all(x -> 0 < x < 1, X)
        @test exactness_error(X, w, 9, :le) < 1e-14
    end

    @testset "normalize = false" begin
        # GH: weight exp(-|x|²); ∫ x₁² x₂⁴ exp(-|x|²) dx = (√π/2)(3√π/4)
        X, w = ghpos(2, 7; normalize = false)
        @test sum(w) ≈ π
        @test sum(w .* X[:, 1] .^ 2 .* X[:, 2] .^ 4) ≈ (sqrt(π) / 2) * (3sqrt(π) / 4)
        # the same integral from the tensor product of FastGaussQuadrature's default rules
        x, v = gausshermite(4)
        @test sum(w .* X[:, 1] .^ 2 .* X[:, 2] .^ 4) ≈ sum(v .* x .^ 2) * sum(v .* x .^ 4)
        # Le: ∫ over [-1,1]³
        X, w = lepos(3, 5; normalize = false)
        @test sum(w) ≈ 8
        @test all(x -> -1 < x < 1, X)
        @test sum(w .* X[:, 1] .^ 2 .* X[:, 3] .^ 2) ≈ 8 / 9
        @test abs(sum(w .* X[:, 1] .* X[:, 2] .^ 2)) < 1e-14
    end

    @testset "even degrees and monotone lookup" begin
        @test ghpos(3, 6) == ghpos(3, 7)
        @test lepos(2, 0) == lepos(2, 1)
        @test nnodes(:le, 2, 10) == nnodes(:le, 2, 11)
    end

    @testset "no rule: error unless pragmatic" begin
        pmax(fam, d) = maximum(r.p for r in available(fam) if r.d == d)
        for (fam, f) in ((:gh, ghpos), (:le, lepos))
            p = pmax(fam, 3) + 2
            @test_throws ArgumentError f(3, p)
            @test_throws ArgumentError nnodes(fam, 3, p)
            @test_throws ArgumentError f(6, 5)                  # no stored rules in six dimensions
            @test nnodes(fam, 3, p; pragmatic = true) ≤ ((p + 1) ÷ 2)^3
        end
        @test_throws ArgumentError ghpos(0, 3)
        @test_throws ArgumentError lepos(2, -1)
        err = try ghpos(6, 5) catch e; e end
        @test occursin("pragmatic = true", err.msg)
    end

    @testset "pragmatic fallback" begin
        # a stored cell is returned unchanged
        @test ghpos(4, 9; pragmatic = true) == ghpos(4, 9)
        @test lepos(2, 21; pragmatic = true) == lepos(2, 21)
        # beyond the stored degrees: a valid rule, cheaper than the product grid
        p = maximum(r.p for r in available(:gh) if r.d == 3) + 2
        X, w = ghpos(3, p; pragmatic = true)
        @test size(X) == (nnodes(:gh, 3, p; pragmatic = true), 3)
        @test all(>(0), w) && sum(w) ≈ 1
        @test exactness_error(X, w, p, :gh) < GATE
        @test size(X, 1) < ((p + 1) ÷ 2)^3
        # beyond the stored dimensions
        for (fam, f) in ((:gh, ghpos), (:le, lepos)), (d, q) in ((6, 7), (7, 5), (8, 3))
            X, w = f(d, q; pragmatic = true)
            n = nnodes(fam, d, q; pragmatic = true)
            @test size(X) == (n, d)
            @test all(>(0), w) && sum(w) ≈ 1
            @test exactness_error(X, w, q, fam) < GATE
            parts = ruleinfo(fam, d, q; pragmatic = true)
            @test sum(r.d for r in parts) == d && prod(r.n for r in parts) == n
            # cheapest: no two-way split of stored rules does better
            @test all(n ≤ nnodes(fam, j, q; pragmatic = true) * nnodes(fam, d - j, q; pragmatic = true) for j in 1:d÷2)
        end
        # normalize = false composes with the fallback
        X, w = lepos(6, 3; normalize = false, pragmatic = true)
        @test sum(w) ≈ 2^6 && sum(w .* X[:, 6] .^ 2) ≈ 2^6 / 3
        # absurd requests fail cleanly instead of exhausting memory
        @test nnodes(:gh, 60, 31; pragmatic = true) isa BigInt
        @test_throws ArgumentError ghpos(60, 31; pragmatic = true)
    end

    @testset "results are copies" begin
        X, w = ghpos(2, 5)
        X .= 0; w .= 0
        X2, w2 = ghpos(2, 5)
        @test sum(w2) ≈ 1 && any(≠(0), X2)
    end

    @testset "extended precision check" begin
        X, w = lepos(2, 9)
        @test exactness_error(big.(X), big.(w), 9, :le) < 1e-15
    end
end
