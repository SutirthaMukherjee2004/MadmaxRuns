@testset "spectrum mode: calibration, check, tau scan, null validation" begin
    root = dirname(@__DIR__)
    isdefined(Main, :validate_main) || Base.include(Main, joinpath(root, "scripts", "validate_nulls.jl"))
    mktempdir() do dir
        # mock spectrum, written as a two-column CSV like a measured one
        cfg = PipelineConfig()
        x, = simulate_spectrum(Xoshiro(11), cfg; model = Heteroscedastic, level = 0.05)
        f = collect(frequencies(cfg.grid))
        csv = joinpath(dir, "spec.csv")
        open(io -> MadmaxRuns.writedlm(io, hcat(f, x), ','), csv, "w")
        out = joinpath(dir, "res.jls")
        Base.invokelatest(Main.main, ["input=$csv", "nbkg=128", "check=true", "taus=1e-7,1e-3",
                                      "output_root=$(joinpath(dir, "out"))", "out=$out"])
        r = load_result(out)
        @test r.result.T > 0 && 0 < r.result.p_pk <= 1
        @test r.calibration.nbkg == 128 && 0 < r.calibration.p_bkg <= 1
        @test r.calibration.check !== nothing && length(r.calibration.check.rows) == 11
        @test r.calibration.check.smooth.eps_max < 0.01
        @test r.calibration.scan !== nothing && length(r.calibration.scan.fp_tau) == 2
        @test 0 < r.calibration.scan.p_scan <= 1
        @test length(r.tauscan) == 2 && all(row -> 0 < row.p <= 1, r.tauscan)
        @test f[1] <= r.run_frequency_hz[1] <= r.run_frequency_hz[2] <= f[end]
        @test r.summary.status isa String
        # unknown option and check without calibration are refused
        @test_throws ArgumentError Base.invokelatest(Main.main, ["input=$csv", "bogus=1", "output_root=$(joinpath(dir, "out"))"])
        @test_throws ArgumentError Base.invokelatest(Main.main, ["input=$csv", "calibrate=false", "check=true",
                                                                 "output_root=$(joinpath(dir, "out"))"])
        # only genuine noise-fit failures are counted, anything else propagates
        @test Main.is_fit_failure(ArgumentError("noise fit did not converge after 60 iterations"))
        @test !Main.is_fit_failure(MethodError(+, ()))
        # null validation on a small directory of signal-free spectra
        nd = joinpath(dir, "nulls"); mkpath(nd)
        for i in 1:6
            xi, = simulate_spectrum(Xoshiro(100 + i), cfg; model = Homoscedastic, level = 0.05)
            open(io -> MadmaxRuns.writedlm(io, hcat(f, xi), ','), joinpath(nd, "s$i.csv"), "w")
        end
        v = Base.invokelatest(Main.validate_main, ["dir=$nd", "output_root=$(joinpath(dir, "out"))"])
        @test length(v.p) == 6 && all(s -> s == "ok", v.status) && all(p -> 0 < p <= 1, v.p)
        @test isfile(joinpath(dir, "out", "data", "null_validation.csv"))
    end
end
