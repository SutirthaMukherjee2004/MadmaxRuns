@testset "boundary term: SVD and removal" begin
    cfg = PipelineConfig(n = 4, m = 64, tau = 1e-6)
    S = filter_setup(cfg; with_pk = false)
    sigma = fill(3.0, S.N)
    a = 5; x = zeros(S.nfull); x[S.J[a]] = 3.0
    c = whiten(linear_residual(S.F, x), S.G, S.keep, sigma)
    col = c[S.rows] ./ sqrt(S.kappa)
    B = boundary_patterns(S.H, fill(3.0, length(S.J)), sigma[S.rows], S.kappa; K = 6)
    Mt = S.H .* (3.0 / sqrt(S.kappa)) ./ 3.0
    @test isapprox(Mt[:, a], col)
    @test isapprox(B.V' * B.V, I(6); atol = 1e-8)
    @test issorted(B.s; rev = true)
    @test B.uncaptured[6] < 1e-5
    xs = 3.0 .* randn(Xoshiro(21), S.nfull)
    z = whiten(linear_residual(S.F, xs), S.G, S.keep, sigma)[S.rows] ./ sqrt(S.kappa)
    zp = remove_boundary(z, B)
    @test squares_T(zp) < squares_T(z)
    @test maximum(abs, B.V' * zp) < 1e-10
end
