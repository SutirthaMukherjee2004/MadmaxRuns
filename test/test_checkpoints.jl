@testset "output folders and live logs" begin
    mktempdir() do dir
        paths = prepare_output_layout(joinpath(dir, "custom results"))
        @test all(isdir, (paths.logs, paths.data, paths.paperplots))
        @test paths == prepare_output_layout(paths.root)
        pipeline_log = joinpath(paths.logs, "pipeline.log")
        inner_log = joinpath(paths.logs, "inner.log")
        with_run_log(pipeline_log) do
            log_message("stage start")
            @test occursin("stage start", read(pipeline_log, String))
            @test_throws ErrorException with_run_log(() -> error("test stage failure"), inner_log)
            @test occursin("Run stopped: test stage failure", read(inner_log, String))
            log_message("stage resumed")
            @test occursin("stage resumed", read(pipeline_log, String))
        end
    end
end

@testset "checkpoints resume only validated stages" begin
    mktempdir() do dir
        checkpoint(x) = MadmaxRuns.CheckpointRun(dir, (; x))
        run = checkpoint(1)
        @test MadmaxRuns.run_stage(() -> [1, 2], run, :input) == [1, 2]
        @test_throws ErrorException MadmaxRuns.run_stage(() -> error("deliberate failure"), run, :noise_fit)
        @test !isfile(joinpath(run.dir, "noise_fit.jls.gz"))

        replay = checkpoint(1)
        @test MadmaxRuns.run_stage(() -> error("input should be reused"), replay, :input) == [1, 2]
        @test MadmaxRuns.run_stage(() -> :fixed, replay, :noise_fit) == :fixed
        @test MadmaxRuns.run_stage(() -> :dependent, replay, :statistics) == :dependent
        @test checkpoint(2).dir != replay.dir

        write(joinpath(replay.dir, "noise_fit.jls.gz"), "corrupt")
        recovery = checkpoint(1)
        @test MadmaxRuns.run_stage(() -> error("input should be reused"), recovery, :input) == [1, 2]
        @test MadmaxRuns.run_stage(() -> :recalculated, recovery, :noise_fit) == :recalculated
        restart = checkpoint(1)
        @test MadmaxRuns.run_stage(() -> error("input should be reused"), restart, :input) == [1, 2]
        @test MadmaxRuns.run_stage(() -> error("fit should be reused"), restart, :noise_fit) == :recalculated
        @test MadmaxRuns.run_stage(() -> :new_statistics, restart, :statistics) == :new_statistics
    end
end

@testset "invalid noise fits stop before scoring" begin
    cfg = PipelineConfig(n = 2, m = 4, W = 8, grid = SpectrumGrid(48, 1e6, 4.218e6), qmax = 12)
    S = filter_setup(cfg; with_pk = false)
    x, = simulate_spectrum(Xoshiro(100), cfg; model = Heteroscedastic, level = 0.05)
    @test_throws ArgumentError run_chain(S, replace(x, x[1] => NaN))
    invalid(sigma, converged) = NoiseFit(1.0, 0.0, zeros(S.N), sigma, 1, converged, false,
                                         [1.0], [0.0], :window, cfg.W)
    @test_throws ArgumentError run_chain(S, x; noise = invalid(ones(S.N), false))
    @test_throws ArgumentError run_chain(S, x; noise = invalid(zeros(S.N), true))
    @test_throws DimensionMismatch run_chain(S, x; noise = invalid(ones(S.N - 1), true))
    mktempdir() do dir
        checkpoint = MadmaxRuns.CheckpointRun(dir, (; cfg, x))
        @test_throws ArgumentError MadmaxRuns._run_chain(S, x, checkpoint;
                                                         noise = invalid(ones(S.N), false))
        @test isfile(joinpath(checkpoint.dir, "noise_fit.failed.jls.gz"))
        @test !isfile(joinpath(checkpoint.dir, "whitening.jls.gz"))
    end
    result = run_chain(S, x)
    @test result.noise.converged && 0 <= result.p_pk <= 1
    @test isnan(run_chain(S, x; compute_pvalues = false).p_pk)
end

module PipelineCLITest
include(joinpath(@__DIR__, "..", "scripts", "run_pipeline.jl"))
end

@testset "measured spectrum CSV validation" begin
    mktempdir() do dir
        csv = joinpath(dir, "spectrum.csv")
        read_csv(text) = (write(csv, text); PipelineCLITest.read_spectrum(csv))
        power, grid = read_csv("4218000,68.12\n4220000,67.94\n4222000,68.31\n")
        @test power == [68.12, 67.94, 68.31]
        @test grid.nbins == 3 && grid.df_hz == 2e3 && grid.fref_hz == 4.218e6
        @test_throws ArgumentError read_csv("1,1\n2,2\n4,3\n")
        @test_throws ArgumentError read_csv("1,1\n1,2\n2,3\n")
        @test_throws ArgumentError read_csv("1,1\n2,NaN\n")
        @test_throws ArgumentError read_csv("frequency_hz,power\n1,2\n")
        @test_throws ArgumentError read_csv("1,2,3\n2,3,4\n")
        @test_throws ArgumentError read_csv("1\n2\n3\n")
        power, grid = read_csv(repeat("35.0\n", PAPER_GRID.nbins))
        @test length(power) == PAPER_GRID.nbins && grid == PAPER_GRID
    end
end
