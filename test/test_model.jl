@testset "step 1: background, noise and line shapes" begin
    b = nominal_background()
    @test length(b) == 25_000
    @test all(isfinite, b)
    bmed = median_background()
    @test 60 < bmed < 80
    for lv in (0.0012, 0.01, 0.05, 0.173)
        s0, s1 = noise_parameters(Homoscedastic, lv, bmed)
        @test s1 == 0 && isapprox(s0, lv * bmed)
        s0, s1 = noise_parameters(Heteroscedastic, lv, bmed)
        @test isapprox(s0, lv * bmed / sqrt(5))
        @test isapprox(s0 + s1 * bmed, lv * bmed)
    end
    for sh in ALL_SHAPES
        p = line_profile(sh, 20.0, 4.0, 4001, 2000.5)
        @test isapprox(sum(p), 20.0; rtol = (sh isa GaussianLine || sh isa TopHatLine ? 1e-6 : 2e-3))
        x = range(-60, 60, length = 240_001); dx = step(x)
        d = diff(line_cdf.(Ref(sh), x, 4.0)) ./ dx
        h = maximum(d) / 2; idx = findall(>=(h), d)
        fw = (idx[end] - idx[1]) * dx
        @test isapprox(fw, 2sqrt(2log(2)) * 4.0; rtol = 2e-3)
    end
    A = line_area_for_peak(AxionLine(), 3.4, 3.0)
    @test isapprox(maximum(line_profile(AxionLine(), A, 3.4, 2001, 1000.5)), 3.0; rtol = 1e-9)
end
