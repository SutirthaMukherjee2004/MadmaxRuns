using MadmaxRuns: parse_cli, paper_tau_grid, autocorrelation, tsum, parabolic_minimum, SpectrumSet, fit_set,
                  eval_set, tau_step, data_generator, case_id, StudyOptions, base_config, filter_setup, with_tau, with

@testset "studies: helpers" begin
    kw = parse_cli(["n=6", "m=48", "levels=0.0012,0.05", "models=homo,hetero", "quick=true", "shapes=axion", "nmc=400", "ms=8:8:32"])
    @test kw[:n] == 6 && kw[:m] == 48 && kw[:levels] == [0.0012, 0.05] && kw[:quick] === true
    @test kw[:models] == [Homoscedastic, Heteroscedastic] && kw[:shapes] == ["axion"]
    @test kw[:nmc] == 400 && kw[:ms] == 8:8:32
    g = paper_tau_grid()
    @test length(g) == 63 && isapprox(g[1], 1e-7) && isapprox(g[end], 0.5) && issorted(g)
    z = randn(Xoshiro(1), 4000)
    rho = autocorrelation(z, 5)
    @test isapprox(rho[1], sum(z[1:end-1] .* z[2:end]) / sum(abs2, z))
    @test tsum(i -> ([i, 2i],), 10)[1] == [55, 110]
    @test isapprox(parabolic_minimum([0.0, 1.0, 2.0], [1.0, 0.0, 1.0]), 10.0)
end

@testset "studies: \u03c3\u0302 reused across \u03c4 equals the full chain" begin
    o = StudyOptions(nspec = 2)
    cfg = with(base_config(o); boundary_removal = true)
    S = filter_setup(with_tau(cfg, 1e-4))
    gen = data_generator(o, cfg, case_id(Heteroscedastic, 0.05); model = Heteroscedastic, level = 0.05)
    set = fit_set(S, "t", gen, 2)
    T, Tp = eval_set(S, set)
    for i in 1:2
        r = run_chain(S, gen(i))
        @test Tp[i] == r.T
    end
end
