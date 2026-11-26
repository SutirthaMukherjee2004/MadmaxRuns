@testset "steps 9\u201310: SQUARES statistic and P_K" begin
    z = randn(Xoshiro(3), 5000)
    T, i0, i1 = squares_T_run(z)
    @test isapprox(T, squares_T(z))
    @test all(z[i0:i1] .>= 0) && isapprox(sum(abs2, z[i0:i1]), T)
    for (t, L) in ((0.1, 1), (3.0, 7), (10.0, 50))
        @test isapprox(pbc_cdf(t, L), MadmaxRuns.RunStatistics.squares_cdf(t, L); atol = 1e-10)
    end
    @test isapprox(pbc_cdf(40.0, 200), pbc_cdf_recursion(40.0, 200); atol = 1e-10)
    @test isapprox(pbc_cdf(25.4, 10_000), MadmaxRuns.pbc_cdf_blocks(25.4, 10_000); atol = 1e-10)
    @test isapprox(pbc_pvalue(25.4, 10_000), 0.115; atol = 0.0005)
    @test isapprox(pbc_cdf(25.4, 24_744), pbc_cdf_recursion(25.4, 24_744); atol = 3e-12, rtol = 0)
    @test 0 < pbc_pvalue(150.0, 24_744) < pbc_pvalue(100.0, 24_744)
    @test_throws DomainError pbc_pvalue(NaN, 100)
    for (h, k, rho) in ((0.3, -0.2, 0.5), (1.2, 0.7, -0.6), (-0.5, 2.0, 0.9))
        ref = let n = 400_000, rng = Xoshiro(1)
            c = 0
            for _ in 1:n
                a = randn(rng); b = rho * a + sqrt(1 - rho^2) * randn(rng)
                c += (a < h) & (b < k)
            end
            c / n
        end
        @test isapprox(Phi2(h, k, rho), ref; atol = 3e-3)
    end
    @test isapprox(Phi2(0.4, -0.1, 0.0), Phi(0.4) * Phi(-0.1))
    hm(k) = (1.0, sqrt(2 / pi), 1.0, 2sqrt(2 / pi), 3.0)[k+1]
    function emom(idx, q)
        s = sort(collect(idx)); p = 1.0; i = 1; K = length(idx)
        while i <= K
            j = i
            while j < K && s[j+1] == s[i]; j += 1; end
            p *= hm(j - i + 1); i = j + 1
        end
        return p / MadmaxRuns.chisq_moment(q, K / 2)
    end
    for q in (1, 3, 5)
        A = Symmetric(randn(Xoshiro(q), q, q)) |> Matrix; b = randn(Xoshiro(10q), q); c = randn(Xoshiro(20q), q)
        m = MadmaxRuns.direction_moments(A, b, c)
        M2 = [emom((i, j), q) for i in 1:q, j in 1:q]; m1 = emom((1,), q)
        Eu = sum(A .* M2); Eb = m1 * sum(b)
        Euu = sum(A[i, j] * A[k, l] * emom((i, j, k, l), q) for i in 1:q, j in 1:q, k in 1:q, l in 1:q)
        Eub = sum(A[i, j] * b[k] * emom((i, j, k), q) for i in 1:q, j in 1:q, k in 1:q)
        @test isapprox(m.Eu, Eu)
        @test isapprox(m.Eb, Eb)
        @test isapprox(m.Vuu, Euu - Eu^2)
        @test isapprox(m.Cub, Eub - Eu * Eb)
        @test isapprox(m.Cbb, b' * M2 * b - Eb^2; atol = 1e-12)
    end
    M0 = MadmaxRuns.PKModel(24872, 24744, Int[]; qmax = 40)
    @test isapprox(pk_cdf(M0, 25.0), pbc_cdf(25.0, 24744))
    @test isapprox(pk_pvalue(M0, 150.0), pbc_pvalue(150.0, 24744); rtol = 1e-12)
    @test_throws DomainError pk_pvalue(M0, 25.0; Eqf = (q, t) -> NaN)
    @test_throws DomainError pk_pvalue(M0, 25.0; Eqf = (q, t) -> 0.0, P0 = 1.0)
    tab = MadmaxRuns.survival_table(t -> t < 80.0 ? pbc_pvalue(t, 10_000) : NaN;
                                    tmax = 80.0, dt = 5.0, allow_invalid = true)
    @test isapprox(pvalue(tab, 25.0), pbc_pvalue(25.0, 10_000))
    @test_throws DomainError pvalue(tab, 80.0)
    finite_tab = PKTable([0.0, 2.0], [0.0, log(0.01)])
    @test isapprox(pvalue(finite_tab, MadmaxRuns.decision_threshold(finite_tab, 0.05)), 0.05)
    @test_throws DomainError pvalue(finite_tab, 3.0)
    @test MadmaxRuns.tabulated_pvalue(_ -> 0.001, finite_tab, 3.0) == 0.001
    K = collect(0:194)
    M1 = MadmaxRuns.PKModel(24872, 24744, K; qmax = 60, order = 1)
    M2m = MadmaxRuns.PKModel(24872, 24744, K; qmax = 60, order = 2)
    @test pk_cdf(M1, 30.0) > pbc_cdf(30.0, 24744)
    @test MadmaxRuns.Eq_second(M2m.eq[5], 25.0) >= MadmaxRuns.Eq_first(M1.eq[5], 25.0)
    Eex, se = Eq_exact(M2m.gamma, 4, 20.0; M = 60_000, rng = Xoshiro(2))
    @test isapprox(MadmaxRuns.Eq_second(M2m.eq[4], 20.0), Eex; atol = max(5se, 2e-3))
end

@testset "P_K over the full T range, fast kernels" begin
    G = transfer_function(SGFilter(4, 64), 24872)
    for tau in (1e-7, 3e-4, 0.05, 0.5)
        M = MadmaxRuns.PKModel(24872, 24744, MadmaxRuns.removed_modes(retained_modes(G, tau)))
        ts = collect(1.0:0.5:400.0)
        p = [pk_pvalue(M, t) for t in ts]
        @test all(v -> 0 < v <= 1, p)                                   # defined everywhere, no error
        @test all(p[i+1] <= p[i] * (1 + 1e-12) for i in 1:length(p)-1)  # monotone in T
        @test isapprox(pk_pvalue(M, M.tail_t - 0.5), -expm1(MadmaxRuns.pk_logcdf(M, M.tail_t - 0.5, nothing, nothing, true)); rtol = 1e-12)
        @test isapprox(pk_pvalue(M, M.tail_t), MadmaxRuns.PK_TAIL_SWITCH; rtol = 1e-6)  # switch at p = 1e-4
        r1 = pk_pvalue(M, 60.0) / pbc_pvalue(60.0, 24744); r2 = pk_pvalue(M, 90.0) / pbc_pvalue(90.0, 24744)
        @test isapprox(r1, r2; rtol = 1e-9)                              # frozen correction factor beyond t*
        @test isfinite(pk_log10pvalue(M, 3000.0)) && pk_log10pvalue(M, 3000.0) < -300
    end
    y = randn(Xoshiro(5), 24872); Y = rfft(y)
    @test isapprox(MadmaxRuns.rfft_cached(y), Y; rtol = 1e-13)
    @test isapprox(MadmaxRuns.irfft_cached(Y, 24872), irfft(Y, 24872); rtol = 1e-13)
    Yc = copy(Y); @test MadmaxRuns.irfft_cached(Yc, 24872) == MadmaxRuns.irfft_cached(Y, 24872) && Yc == Y  # input preserved
    B = randn(Xoshiro(6), 7, 5); A = B' * B + I; b = randn(Xoshiro(7), 5)
    @test isapprox(MadmaxRuns.spd_solve!(copy(A), copy(b)), A \ b; rtol = 1e-12)
    @test_throws PosDefException MadmaxRuns.spd_solve!([1.0 2.0; 2.0 1.0], [1.0, 1.0])
    for i in 1:200
        z = randn(Xoshiro(100 + i), 3000)
        @test squares_T_run(z)[1] == squares_T(z)                       # same T as RunStatistics.t_obs
    end
end
