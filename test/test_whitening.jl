@testset "steps 4\u20138: Fourier whitening" begin
    F = SGFilter(4, 64); N = 25_000 - 128
    G = transfer_function(F, N)
    keep = retained_modes(G, 1e-7)
    @test !keep[1] && !keep[end]
    @test MadmaxRuns.n_removed_rfft(keep) == 97
    kappa = retained_fraction(keep, N)
    @test 0.992 < kappa < 0.9925
    @test sum(retained_modes(G, 1e-7; drop_nyquist = false)) == sum(keep) + 1
    @test issorted(.!keep[1:200]; rev = true)
    rng = Xoshiro(5); sigma = 2.0
    v = zeros(N); ns = 300; acc1 = 0.0
    for _ in 1:ns
        n = sigma .* randn(rng, N)
        r = irfft(rfft(n) .* G, N)
        c = whiten(r, G, keep, fill(sigma, N))
        v .+= c .^ 2
        acc1 += mean(c[1:end-1] .* c[2:end])
    end
    @test isapprox(mean(v) / ns, kappa; rtol = 0.01)
    Q1 = -sum(MadmaxRuns.mode_weight(f, N) * cos(2 * pi * f / N) for f in removed_modes(keep)) / N
    @test isapprox(acc1 / ns, Q1; atol = 1e-3)
end
