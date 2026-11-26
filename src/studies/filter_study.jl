function study_filter(o::StudyOptions)
    N = PAPER_GRID.nbins - 2 * 64
    specs = [(2, 32, :SG), (4, 32, :SG), (4, 64, :SG), (4, 32, :SGW)]
    out = NamedTuple[]
    for (n, m, k) in specs
        F = SGFilter(n, m; kind = k)
        G = transfer_function(F, N)
        push!(out, (; n, m, kind = k, d = collect(-m:m), S = F.S, f = collect(0:div(N, 2)) ./ N, G2 = abs2.(G)))
    end
    res = (; filters = out, N)
    save_result(result_path(o, "filter"), res)
    return res
end
