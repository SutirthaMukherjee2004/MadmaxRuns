const TCUT_NTAU = [195, 324, 387]

function study_tcut(o::StudyOptions; ntau = TCUT_NTAU, nmc::Integer = 500_000, Mexact::Integer = 20_000,
                    tgrid = collect(1.0:0.5:40.0), qexact::Integer = 60)
    N = PAPER_GRID.nbins - 2o.m; rows = tested_bins(N, o.m); Ng = length(rows)
    Pbc = [pbc_cdf(t, Ng) for t in tgrid]
    alphas = [0.05, 0.01, 1e-3, 1e-4]
    out = Dict{Int,Any}()
    for nt in ntau
        K = collect(0:nt-1)
        keep = trues(div(N, 2) + 1); keep[K .+ 1] .= false
        kappa = retained_fraction(keep, N)
        T = zeros(nmc)
        Threads.@threads for i in 1:nmc
            Y = rfft(randn(stream_rng(o.seed + 7000 + nt, i), N))
            Y[.!keep] .= 0
            T[i] = squares_T(view(irfft(Y, N), rows) ./ sqrt(kappa))
        end
        sort!(T)
        M1 = PKModel(N, Ng, K; order = 1); M2 = PKModel(N, Ng, K; order = 2)
        qmax = min(qexact, length(M2.eq))
        Etab = zeros(qmax, length(tgrid))
        Threads.@threads for idx in 1:qmax*length(tgrid)
            q, k = (idx - 1) % qmax + 1, div(idx - 1, qmax) + 1
            Etab[q, k] = Eq_exact(M2.gamma, q, tgrid[k]; M = Mexact, rng = stream_rng(o.seed + 9000 + nt, idx))[1]
        end
        P1 = [pk_cdf(M1, t; P0 = Pbc[k]) for (k, t) in enumerate(tgrid)]
        P2 = [pk_cdf(M2, t; P0 = Pbc[k]) for (k, t) in enumerate(tgrid)]
        Pe = [pk_cdf(M2, t; P0 = Pbc[k], Eqf = (q, _) -> q <= qmax ? Etab[q, k] : Eq_second(M2.eq[q], t)) for (k, t) in enumerate(tgrid)]
        Smc = [(length(T) - searchsortedfirst(T, t) + 1) / length(T) for t in tgrid]
        err = sqrt.(Smc .* (1 .- Smc) ./ length(T))
        edges = collect(0.0:0.5:maximum(T) + 0.5)
        dens = [count(x -> edges[j] <= x < edges[j+1], T) / (length(T) * 0.5) for j in 1:length(edges)-1]
        cal = Dict{String,Vector{Float64}}()
        for (name, cdf) in (("exact", nothing), ("second", t -> pk_cdf(M2, t)), ("first", t -> pk_cdf(M1, t)), ("pbc", t -> pbc_cdf(t, Ng)))
            cal[name] = map(alphas) do alpha
                talpha = cdf === nothing ? interp_threshold(tgrid, Pe, alpha) : pk_threshold(cdf, alpha; lo = 0.5, hi = 80.0, tol = 1e-4)
                (length(T) - searchsortedfirst(T, talpha) + 1) / length(T)
            end
        end
        out[nt] = (; K, kappa, T_sample = T[1:max(1, div(length(T), 50000)):end], density = dens, edges,
                   S_mc = Smc, S_err = err, S_exact = 1 .- Pe, S_second = 1 .- P2, S_first = 1 .- P1, S_pbc = 1 .- Pbc,
                   ratio_exact = (1 .- Pe) ./ Smc, ratio_second = (1 .- P2) ./ Smc, ratio_first = (1 .- P1) ./ Smc,
                   band = err ./ Smc, calibration = cal, Etab)
        log_message(@sprintf("tcut n_\u03c4 = %d: S_mc(25) = %.4f, exact %.4f, second %.4f, first %.4f, P_BC %.4f", nt,
                             Smc[argmin(abs.(tgrid .- 25))], 1 - Pe[argmin(abs.(tgrid .- 25))],
                             1 - P2[argmin(abs.(tgrid .- 25))], 1 - P1[argmin(abs.(tgrid .- 25))], 1 - Pbc[argmin(abs.(tgrid .- 25))]))
    end
    res = (; tgrid, N, Ng, cases = out, alphas = alphas, nmc, ntau = collect(ntau))
    save_result(result_path(o, "tcut"), res)
    return res
end

function interp_threshold(t::AbstractVector, P::AbstractVector, alpha::Real)
    S = 1 .- P
    k = findfirst(<=(alpha), S)
    k === nothing && return t[end]
    k == 1 && return t[1]
    a = (log(S[k-1]) - log(alpha)) / (log(S[k-1]) - log(S[k]))
    return t[k-1] + a * (t[k] - t[k-1])
end
