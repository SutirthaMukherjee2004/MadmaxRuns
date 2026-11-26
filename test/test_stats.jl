@testset "steps 9\u201310: SQUARES statistic and P_K" begin
    z = randn(Xoshiro(3), 5000)
    T, i0, i1 = squares_T_run(z)
    @test isapprox(T, squares_T(z))
    @test all(z[i0:i1] .>= 0) && isapprox(sum(abs2, z[i0:i1]), T)
    for (t, L) in ((10.0, 50), (25.0, 150), (40.0, 200))
        @test isapprox(pbc_cdf(t, L), pbc_cdf_recursion(t, L); atol = 1e-10)
    end
    @test isapprox(pbc_cdf(25.4, 10_000), MadmaxRuns.pbc_cdf_blocks(25.4, 10_000); atol = 1e-10)
    @test 0.11 < pbc_pvalue(25.4, 10_000) < 0.12
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
    K = collect(0:194)
    M1 = MadmaxRuns.PKModel(24872, 24744, K; qmax = 60, order = 1)
    M2m = MadmaxRuns.PKModel(24872, 24744, K; qmax = 60, order = 2)
    @test pk_cdf(M1, 30.0) > pbc_cdf(30.0, 24744)
    @test MadmaxRuns.Eq_second(M2m.eq[5], 25.0) >= MadmaxRuns.Eq_first(M1.eq[5], 25.0)
    Eex, se = Eq_exact(M2m.gamma, 4, 20.0; M = 60_000, rng = Xoshiro(2))
    @test isapprox(MadmaxRuns.Eq_second(M2m.eq[4], 20.0), Eex; atol = max(5se, 2e-3))
end
