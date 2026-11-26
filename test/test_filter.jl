@testset "step 2: Savitzky\u2013Golay filter" begin
    for (n, m) in ((2, 32), (4, 32), (4, 64), (6, 64))
        F = SGFilter(n, m)
        @test isapprox(sum(F.S), 1; atol = 1e-12)
        @test isapprox(noise_pass_fraction(F), 1 - F.S[m+1]; atol = 1e-12)
        for p in 1:n
            @test abs(kernel_moments(F, p)) < 1e-6 * m^p
        end
        pst = iseven(n) ? n + 2 : n + 1
        @test abs(kernel_moments(F, pst)) > 1e-3 * m^pst
        G = transfer_function(F, 4096)
        @test abs(G[1]) < 1e-12
        x = [(i / 500)^n - 0.3(i / 500) + 2 for i in 1:3000]
        @test maximum(abs, linear_residual(F, x)) < 1e-8
        y = randn(Xoshiro(1), 1000)
        r = linear_residual(F, y)
        @test isapprox(r[10], y[10+m] - sum(F.S[k+m+1] * y[10+m+k] for k in -m:m))
        N = 1024; per = randn(Xoshiro(2), N)
        xp = [per[mod1(i - m, N)] for i in 1:N+2m]
        rc = linear_residual(F, xp)
        Gc = transfer_function(F, N)
        @test maximum(abs, rc .- irfft(rfft(per) .* Gc, N)) < 1e-10
    end
    F = SGFilter(4, 32; kind = :SGW)
    @test isapprox(sum(F.S), 1; atol = 1e-12)
    @testset "weighted SG: FFT form = direct form" begin
        F = SGFilter(4, 16); nf = 600
        x = randn(Xoshiro(3), nf) .+ 50; wt = 1 ./ (1 .+ 0.3rand(Xoshiro(4), nf)) .^ 2
        c = MadmaxRuns.WSGCache(F, nf)
        @test maximum(abs, weighted_smooth(c, x, wt) .- weighted_smooth_direct(F, x, wt)) < 1e-8
        @test maximum(abs, weighted_smooth(c, x, ones(nf)) .- smooth_unweighted(F, x)) < 1e-8
    end
end
