@testset "step 3: iterative noise estimate" begin
    cfg = PipelineConfig(n = 4, m = 64)
    F = SGFilter(4, 64); c = MadmaxRuns.WSGCache(F, 25_000)
    bmed = median_background()
    for (model, lv) in ((Heteroscedastic, 0.05), (Homoscedastic, 0.05), (Heteroscedastic, 0.0012))
        x, b, sigma, _ = simulate_spectrum(Xoshiro(7), cfg; model, level = lv, bmed)
        s0, s1 = noise_parameters(model, lv, bmed)
        for est in (:window, :lsq)
            fit = estimate_noise(F, c, x; W = 96, estimator = est)
            @test fit.converged
            @test fit.iterations <= 20
            @test length(fit.sigma0_hist) == fit.iterations + 1
            sigmat = sigma[65:end-64]
            @test median(abs.(fit.sigma ./ sigmat .- 1)) < 0.03
            model == Homoscedastic && @test abs(fit.sigma1) * bmed < 0.03 * s0
            model == Heteroscedastic && @test isapprox(fit.sigma1, s1; rtol = 0.08)
        end
    end
    rng = Xoshiro(9); bl = collect(range(10, 100, length = 96 * 300))
    sigmal = 0.5 .+ 0.02 .* bl
    s0, s1, _, _ = window_fit(bl .+ sigmal .* randn(rng, length(bl)), 96, ones(length(bl)))
    @test isapprox(s0, 0.5; rtol = 0.1)
    @test isapprox(s1, 0.02; rtol = 0.05)
end
