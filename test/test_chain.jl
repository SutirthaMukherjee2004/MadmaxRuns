@testset "full chain, steps 1\u201310, and the null hypotheses" begin
    cfg = PipelineConfig(n = 4, m = 64, tau = 1e-7)
    S = filter_setup(cfg)
    bmed = median_background()
    @test S.Ng == 24_744
    A = line_area_for_peak(AxionLine(), 3.4, 3.0)
    for model in (Homoscedastic, Heteroscedastic)
        x, b, sigma, s = simulate_spectrum(Xoshiro(100), cfg; model, level = 0.05, bmed,
                                       signal = (AxionLine(), A, 3.4, 10.5e6))
        res = run_chain(S, x)
        @test res.noise.converged
        @test length(res.c) == S.N && length(res.zp) == S.Ng
        @test res.p_pk < 1e-3
        lo = frequency_bin(cfg.grid, 10.5e6) - 64
        @test lo - 6 <= res.run[1] + 64 <= lo + 6
        x0, _, _, _ = simulate_spectrum(Xoshiro(101), cfg; model, level = 0.05, bmed)
        r0 = run_chain(S, x0)
        @test abs(mean(r0.zp)) < 0.05 && 0.9 < var(r0.zp) < 1.1
    end
    iid = null_iid(S, 64; seed = 1)
    @test issorted(iid.T) && length(iid.T) == 64
    @test 15 < median(iid.T) < 35
    x0, _, sigma0, _ = simulate_spectrum(Xoshiro(5), cfg; model = Heteroscedastic, level = 0.05, bmed)
    r0 = run_chain(S, x0)
    filt = null_filt(S, [r0.noise.sigma], 8; seed = 2)
    @test length(filt.T) == 8 && all(isfinite, filt.T)
    @test 0 < mc_pvalue(r0.T, iid) <= 1
end
