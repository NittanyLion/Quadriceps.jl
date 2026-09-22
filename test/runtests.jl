using Quadriceps, Test
using Aqua
using FastGaussQuadrature: gausshermite, gausslegendre
using Quadriceps: available, nnodes, ruleinfo, exactness_error

# Every stored rule passed a 1e-11 gate when the data were built; the catalog records the error.
const GATE = 1e-11

@testset "Quadriceps" begin

    @testset "Aqua" begin
        Aqua.test_all(Quadriceps)
    end

    @testset "catalog" begin
        @test !isempty(available(:gh)) && !isempty(available(:le))
        @test all(r -> isodd(r.p) && r.d ≥ 2 && r.n ≥ 1, vcat(available(:gh), available(:le)))
        @test_throws ArgumentError available(:la)
    end

    @testset "stored $fam rules" for (fam, f) in ((:gh, ghpos), (:le, lepos))
        for r in available(fam)
            X, w = f(r.d; p = r.p)
            @test f(r.d, (r.p + 1) ÷ 2) == (X, w)               # q form: p = 2q - 1
            @test size(X) == (r.n, r.d) && length(w) == r.n
            @test all(>(0), w)
            @test sum(w) ≈ 1
            @test exactness_error(X, w, r.p, fam) < GATE
            fam ≡ :le && @test all(x -> 0 < x < 1, X)
            @test nnodes(fam, r.d; p = r.p) ≤ r.n               # a higher degree may be cheaper, never dearer
            @test nnodes(fam, r.d; p = r.p, pragmatic = true) == nnodes(fam, r.d; p = r.p)   # no product beats a stored rule
        end
    end

    @testset "one dimension" begin
        x, v = gausshermite(4; normalize = true)
        X, w = ghpos(1, 4)
        @test size(X) == (4, 1) && vec(X) == x && w == v
        x, v = gausshermite(4)
        X, w = ghpos(1; p = 6, normalize = false)           # even degree: next odd
        @test vec(X) ≈ x && w ≈ v
        x, v = gausslegendre(5)
        X, w = lepos(1, 5; normalize = false)
        @test vec(X) ≈ x && w ≈ v
        X, w = lepos(1, 5)
        @test sum(w) ≈ 1 && all(x -> 0 < x < 1, X)
        @test exactness_error(X, w, 9, :le) < 1e-14
    end

    @testset "normalize = false" begin
        # GH: weight exp(-|x|²); ∫ x₁² x₂⁴ exp(-|x|²) dx = (√π/2)(3√π/4)
        X, w = ghpos(2, 4; normalize = false)
        @test sum(w) ≈ π
        @test sum(w .* X[:, 1] .^ 2 .* X[:, 2] .^ 4) ≈ (sqrt(π) / 2) * (3sqrt(π) / 4)
        # the same integral from the tensor product of FastGaussQuadrature's default rules
        x, v = gausshermite(4)
        @test sum(w .* X[:, 1] .^ 2 .* X[:, 2] .^ 4) ≈ sum(v .* x .^ 2) * sum(v .* x .^ 4)
        # Le: ∫ over [-1,1]³
        X, w = lepos(3, 3; normalize = false)
        @test sum(w) ≈ 8
        @test all(x -> -1 < x < 1, X)
        @test sum(w .* X[:, 1] .^ 2 .* X[:, 3] .^ 2) ≈ 8 / 9
        @test abs(sum(w .* X[:, 1] .* X[:, 2] .^ 2)) < 1e-14
    end

    @testset "even degrees and monotone lookup" begin
        @test ghpos(3; p = 6) == ghpos(3; p = 7) == ghpos(3, 4)
        @test lepos(2; p = 0) == lepos(2, 1)
        @test nnodes(:le, 2; p = 10) == nnodes(:le, 2, 6)
    end

    @testset "no rule: error unless pragmatic" begin
        qmax(fam, d) = maximum((r.p + 1) ÷ 2 for r in available(fam) if r.d == d)
        for (fam, f) in ((:gh, ghpos), (:le, lepos))
            q = qmax(fam, 3) + 1
            @test_throws ArgumentError f(3, q)
            @test_throws ArgumentError f(3; p = 2q - 2)
            @test_throws ArgumentError nnodes(fam, 3, q)
            @test_throws ArgumentError f(6, 3)                  # no stored rules in six dimensions
            @test nnodes(fam, 3, q; pragmatic = true) ≤ q^3
        end
        @test_throws ArgumentError ghpos(0, 3)
        @test_throws ArgumentError lepos(2, 0)
        @test_throws ArgumentError lepos(2; p = -1)
        @test_throws UndefKeywordError ghpos(2)
        err = try ghpos(6, 3) catch e; e end
        @test occursin("pragmatic = true", err.msg)
    end

    @testset "pragmatic fallback" begin
        # a stored cell is returned unchanged
        @test ghpos(4, 5; pragmatic = true) == ghpos(4, 5)
        @test lepos(2, 11; pragmatic = true) == lepos(2, 11)
        # beyond the stored degrees: a valid rule, cheaper than the product grid (d = 4: the split 2 + 2
        # into stored d = 2 rules beats the grid; at d = 3 the stored d = 2 rules end at the same degree as
        # the d = 3 ones, so the fallback there is the bare grid)
        q = maximum((r.p + 1) ÷ 2 for r in available(:gh) if r.d == 4) + 1; p = 2q - 1
        X, w = ghpos(4, q; pragmatic = true)
        @test size(X) == (nnodes(:gh, 4, q; pragmatic = true), 4)
        @test all(>(0), w) && sum(w) ≈ 1
        @test exactness_error(X, w, p, :gh) < GATE
        @test size(X, 1) < q^4
        @test nnodes(:gh, 3, 18; pragmatic = true) == 18^3                # d = 3, q = 18: nothing beats the grid
        # beyond the stored dimensions
        for (fam, f) in ((:gh, ghpos), (:le, lepos)), (d, q) in ((6, 4), (7, 3), (8, 2))
            X, w = f(d, q; pragmatic = true)
            n = nnodes(fam, d, q; pragmatic = true)
            @test size(X) == (n, d)
            @test all(>(0), w) && sum(w) ≈ 1
            @test exactness_error(X, w, 2q - 1, fam) < GATE
            parts = ruleinfo(fam, d, q; pragmatic = true)
            @test sum(r.d for r in parts) == d && prod(r.n for r in parts) == n
            # cheapest: no two-way split of stored rules does better
            @test all(n ≤ nnodes(fam, j, q; pragmatic = true) * nnodes(fam, d - j, q; pragmatic = true) for j in 1:d÷2)
        end
        # normalize = false composes with the fallback
        X, w = lepos(6, 2; normalize = false, pragmatic = true)
        @test sum(w) ≈ 2^6 && sum(w .* X[:, 6] .^ 2) ≈ 2^6 / 3
        # absurd requests fail cleanly instead of exhausting memory
        @test nnodes(:gh, 60, 16; pragmatic = true) isa BigInt
        @test_throws ArgumentError ghpos(60, 16; pragmatic = true)
    end

    @testset "data file" begin
        using Quadriceps: readbinindex, readblock, writebin, binpath, BIN, INDEX
        @test startswith(readline(binpath()), "QUADRICEPS1 fmt=1 endian=little cells=$(length(INDEX)) index_fields=8 float=binary64")
        @test keys(BIN) == keys(INDEX)
        @test all(BIN[k].n == r.n && BIN[k].source_id == r.source_id && BIN[k].nbytes == r.n * (r.d + 1) * 8 for (k, r) in INDEX)
        @test filesize(binpath()) == maximum(e.offset + e.nbytes for e in values(BIN))      # no slack
        @test isempty([f for (_, _, fs) in walkdir(dirname(binpath())) for f in fs if endswith(f, ".csv")])
        # round trip through a second file
        rules = Dict(k => (ghpos(k[2]; p = k[3])..., 7) for k in [(:gh, 2, 5), (:gh, 3, 7)])
        rules[(:le, 2, 9)] = (lepos(2; p = 9)..., 16)
        mktempdir() do dir
            path = writebin(joinpath(dir, "t.bin"), rules)
            bin = readbinindex(path)
            @test keys(bin) == keys(rules)
            @test all(readblock(path, bin[k], k[2]) == rules[k][1:2] && bin[k].source_id == rules[k][3] for k in keys(rules))
        end
    end

    @testset "results are copies" begin
        X, w = ghpos(2, 3)
        X .= 0; w .= 0
        X2, w2 = ghpos(2, 3)
        @test sum(w2) ≈ 1 && any(≠(0), X2)
    end

    @testset "extended precision check" begin
        X, w = lepos(2, 5)
        @test exactness_error(big.(X), big.(w), 9, :le) < 1e-15
    end

    @testset "binary128 codec" begin
        enc(x) = Quadriceps.encode128(BigFloat(x; precision = 113))
        @test enc(1) == UInt128(0x3fff) << 112
        @test enc(-2.5) == (UInt128(1) << 127) | (UInt128(0x4000) << 112) | (UInt128(1) << 110)
        @test enc("0.1") == 0x3ffb999999999999999999999999999a          # the binary128 nearest to 1/10
        @test enc(0) == 0
        for x ∈ (BigFloat(π; precision = 113), -BigFloat("7.0785769767068e-33"; precision = 113), BigFloat(2; precision = 113)^-40)
            @test Quadriceps.decode128(Quadriceps.encode128(x)) == x
        end
        @test_throws ArgumentError Quadriceps.encode128(BigFloat(π; precision = 256))
    end

    @testset "number types beyond Float64: the binary128 data" begin
        eps128 = 2.0^-112
        @test Quadriceps.bits(Float64) == 53 && !Quadriceps.wide(Float64)
        setprecision(BigFloat, 113) do                # ≤ 113 bits: the binary128 data, rules128.bin
            @test !Quadriceps.wide(BigFloat)
            for fam ∈ (:gh, :le)
                ext = Quadriceps.extended(fam)
                @test length(ext) == length(available(fam))
                pos = fam ≡ :gh ? ghpos : lepos
                for r ∈ ext
                    @test r.relerr128 ≤ 10eps128 && r.relerr80 < 1e-66
                    r.n ≤ 300 || continue
                    X, w = pos(BigFloat, r.d; p = r.p)
                    X64, w64 = pos(r.d; p = r.p)
                    @test size(X) == size(X64) && eltype(X) ≡ BigFloat && precision(w[1]) == 113 && all(>(0), w)
                    # the binary128 rule rounds to the Float64 rule: bit for bit, or one unit in the last place off
                    @test all(abs(Float64(a) - b) ≤ eps(b) for (a, b) ∈ zip(X, X64)) && all(abs(Float64(a) - b) ≤ eps(b) for (a, b) ∈ zip(w, w64))
                    err = setprecision(() -> exactness_error(BigFloat.(X), BigFloat.(w), r.p, fam), BigFloat, 320)
                    @test err ≤ 10eps128
                    @test isapprox(Float64(err), r.relerr128; rtol = 1e-2) || r.relerr128 == 0
                end
            end
            # ZENODO_ONLY: GH d=2 stops at p=33 (the deposit's ceiling); beyond it there is nothing stored,
            # and `pragmatic` rebuilds the degree from the 20² / 21² product grid, in any type
            @test_throws ArgumentError ghpos(BigFloat, 2; p = 39)
            @test size(ghpos(BigFloat, 2; p = 39, pragmatic = true)[1], 1) == 400
            @test size(ghpos(BigFloat, 2; p = 41, pragmatic = true)[1], 1) == 441
        end
        @test eltype(ghpos(Float32, 3, 4)[1]) ≡ Float32 && Float32.(ghpos(3, 4)[2]) == ghpos(Float32, 3, 4)[2]
        @test ghpos(Float64, 3, 4) == ghpos(3, 4) && lepos(Float64, 2; p = 9) == lepos(2, 5)
        # every catalog cell is in the binary128 catalog too, so no stored cell is Float64 only
        @test Set((r.d, r.p) for fam ∈ (:gh, :le) for r ∈ Quadriceps.available(fam)) ==
              Set((r.d, r.p) for fam ∈ (:gh, :le) for r ∈ Quadriceps.extended(fam))
        # one-dimensional Gauss rules, refined beyond Float64 (default BigFloat: the wide path)
        for (pos, fam, fgq) ∈ ((ghpos, :gh, q -> gausshermite(q; normalize = true)), (lepos, :le, q -> (g = gausslegendre(q); ((g[1] .+ 1) ./ 2, g[2] ./ 2))))
            for q ∈ (1, 2, 7, 30)
                X, w = pos(BigFloat, 1, q)
                @test exactness_error(X, w, 2q - 1, fam) < 1e-60
                x0, w0 = fgq(q)
                @test Float64.(vec(X)) ≈ x0 && Float64.(w) ≈ w0
            end
        end
        # frames and tensor products in T
        X, w = ghpos(BigFloat, 3, 4; normalize = false)
        @test abs(sum(w) - sqrt(big(π))^3) < 1e-32 && Float64.(X) ≈ ghpos(3, 4; normalize = false)[1]
        X, w = lepos(BigFloat, 2, 5; normalize = false)
        @test abs(sum(w) - 4) < 1e-32 && all(x -> -1 < x < 1, X)
        X, w = ghpos(BigFloat, 7, 3; pragmatic = true)
        @test eltype(w) ≡ BigFloat && size(X) == (nnodes(:gh, 7, 3; pragmatic = true), 7) && abs(sum(w) - 1) < 1e-32
        @test exactness_error(X, w, 5, :gh) < 10eps128
    end

    @testset "binary320 codec" begin
        using Quadriceps: encode320, decode320, PREC80
        enc(x) = encode320(BigFloat(x; precision = PREC80))
        @test enc(1) == (0, 0, 0, 0, UInt64(0x3fff) << 48)                    # exponent field in the top word
        @test enc(0) == (0, 0, 0, 0, 0)
        @test enc(-2.5)[5] == (UInt64(1) << 63) | (UInt64(0x4000) << 48) | (UInt64(1) << 46)   # 1.01b: second fraction bit
        for x ∈ (BigFloat(π; precision = PREC80), -BigFloat("7.0785769767068e-33"; precision = PREC80),
                 BigFloat(2; precision = PREC80)^-40, BigFloat("0." * "1234567890"^8; precision = PREC80))
            @test decode320(enc(x)) == x
        end
        @test_throws ArgumentError encode320(BigFloat(π; precision = 400))
        # the deposit's 80-digit decimals survive the round trip to within the format's last bit
        s = "9.3173579904737320568086554096983364255055931636511679768181393607476690421431501e-01"
        @test abs(decode320(encode320(BigFloat(s; precision = PREC80))) - BigFloat(s; precision = 512)) < BigFloat(2)^-304
    end

    @testset "80 digits: rules80.bin" begin
        using Quadriceps: bin80path, BIN80, BIN128
        @test precision(BigFloat) == 256 && Quadriceps.wide(BigFloat)      # the default BigFloat is wide
        @test startswith(readline(bin80path()), "QUADRICEPS1 fmt=1 endian=little cells=$(length(INDEX)) index_fields=8 float=binary320")
        @test keys(BIN80) == keys(BIN128) == keys(INDEX)
        @test all(BIN80[k].n == e.n && BIN80[k].nbytes == e.n * (k[2] + 1) * 40 for (k, e) ∈ BIN128)
        @test filesize(bin80path()) == maximum(e.offset + e.nbytes for e ∈ values(BIN80))
        for fam ∈ (:gh, :le)
            for r ∈ Quadriceps.extended(fam)
                info = Quadriceps.INDEX[(fam, r.d, r.p)]
                X80, w80 = Quadriceps.stored80(info)
                @test size(X80) == (r.n, r.d) && precision(w80[1]) == Quadriceps.PREC80 && all(>(0), w80)
                @test abs(sum(w80) - 1) < 1e-66                                    # the files carry 80 digits: ~1e-70
                # the 80-digit rule rounds to the binary128 rule exactly (and so to the Float64 rule)
                X128, w128 = Quadriceps.stored128(info)
                @test all(BigFloat(a; precision = 113) == b for (a, b) ∈ zip(X80, X128)) &&
                      all(BigFloat(a; precision = 113) == b for (a, b) ∈ zip(w80, w128))
                fam ≡ :le && @test all(x -> 0 < x < 1, X80)
                if r.n ≤ 300
                    err = setprecision(() -> exactness_error(X80, w80, r.p, fam), BigFloat, 400)
                    @test err < 1e-66
                    @test isapprox(Float64(err), r.relerr80; rtol = 0.2) || err < 1e-75
                end
            end
        end
        # through the API: a default BigFloat is the 80-digit rule rounded to 256 bits …
        X, w = ghpos(BigFloat, 3, 4)
        X80, w80 = Quadriceps.stored80(Quadriceps.INDEX[(:gh, 3, 7)])
        @test precision(w[1]) == 256 && X == BigFloat.(X80) && w == BigFloat.(w80)
        @test setprecision(() -> exactness_error(X, w, 7, :gh), BigFloat, 400) < 1e-66
        # … at 1000 bits it is the same 80-digit rule, not a better one …
        X2, w2 = setprecision(() -> ghpos(BigFloat, 3, 4), BigFloat, 1000)
        @test precision(w2[1]) == 1000 && [BigFloat(x; precision = 256) for x ∈ w2] == w
        # … and at 113 bits it is the package's binary128 rule
        @test setprecision(() -> ghpos(BigFloat, 3, 4)[2], BigFloat, 113) == Quadriceps.stored128(Quadriceps.INDEX[(:gh, 3, 7)])[2]
        # frames and one-dimensional factors keep up with the data
        X, w = lepos(BigFloat, 2, 5; normalize = false)
        @test abs(sum(w) - 4) < 1e-66 && all(x -> -1 < x < 1, X)
        X, w = lepos(BigFloat, 6, 3; pragmatic = true)                      # a product of stored and Gauss factors
        @test size(X, 2) == 6 && abs(sum(w) - 1) < 1e-66 && exactness_error(X, w, 5, :le) < 1e-66
        # results are copies of the cache
        X, w = ghpos(BigFloat, 2, 3); w .= 0
        @test all(>(0), ghpos(BigFloat, 2, 3)[2])
    end
end
