function study_table_a(o::StudyOptions)
    L = NOBKG_L; cen = L / 2
    SH = [parse_shape(s) for s in o.shapes]; AS = o.amplitudes; WS = o.widths
    cutoff = decision_threshold(pbc_table(L; tmax = 150.0), 0.05)
    Tn = zeros(o.nnull)
    Threads.@threads for s in 1:o.nnull
        Tn[s] = squares_T(randn(spectrum_rng(o, 2, s), L))
    end
    null = NullSet(:iid, sort!(Tn))
    Ff = zeros(length(SH), length(AS), length(WS)); Fm = similar(Ff)
    for (i, sh) in enumerate(SH), (j, w) in enumerate(WS)
        prof = line_profile(sh, 1.0, w, L, cen)
        for (a, A) in enumerate(AS)
            T = zeros(o.nspec)
            Threads.@threads for s in 1:o.nspec
                T[s] = squares_T(randn(spectrum_rng(o, 1, s), L) .+ A .* prof)
            end
            Ff[i, a, j] = mean(T .>= cutoff)
            Fm[i, a, j] = mean(mc_pvalue.(T, Ref(null)) .<= 0.05)
        end
        log_message("table (a): $(shape_name(sh)), w = $w: " * join([@sprintf("%.3f", Ff[i, a, j]) for a in eachindex(AS)], " "))
    end
    res = (; part = :a, shapes = shape_name.(SH), AS, WS, formula = Ff, mc = Fm, nspec = o.nspec, nnull = o.nnull)
    save_result(result_path(o, "table_a" * shape_suffix(o)), res)
    write_table_csv(joinpath(o.outdir, "table_a" * shape_suffix(o) * ".csv"), res)
    return res
end

shape_suffix(o::StudyOptions) = length(o.shapes) == 4 && length(o.widths) == 6 ? "" :
    "_" * join(replace.(o.shapes, "-" => ""), "-") * "_w" * join(string.(Int.(round.(o.widths))), "-")

function study_table_b(o::StudyOptions)
    cfg = base_config(o; tau = o.tau, boundary_removal = true)
    S = filter_setup(cfg)
    cutoff = decision_threshold(S.pktab, 0.05)
    bmed = median_background(); model = Homoscedastic; lv = o.level
    SH = [parse_shape(s) for s in o.shapes]; AS = o.amplitudes; WS = o.widths
    log_message(@sprintf("table (b): SG n = %d, m = %d, \u03c4 = %.1e, %d rfft modes removed, \u03ba = %.5f",
                         o.n, o.m, o.tau, n_removed_rfft(S.keep), S.kappa))
    null = null_iid(S, o.nnull; seed = o.seed + 11)
    sigma0, sigma1 = noise_parameters(model, lv, bmed)
    cen = S.nfull / 2
    base = Vector{Vector{Float64}}(undef, o.nspec); sig = similar(base)
    Threads.@threads for s in 1:o.nspec
        rng = spectrum_rng(o, 3, s)
        b = background(rng); sigma = noise_sigma(sigma0, sigma1, b)
        base[s] = b .+ sigma .* randn(rng, S.nfull); sig[s] = sigma
    end
    Ff = zeros(length(SH), length(AS), length(WS)); Fm = similar(Ff)
    T0 = zeros(o.nspec)
    Threads.@threads for s in 1:o.nspec
        T0[s] = run_chain(S, base[s]; compute_pvalues = false).T
    end
    for (i, sh) in enumerate(SH), (j, w) in enumerate(WS)
        prof = line_profile(sh, 1.0, w, S.nfull, cen)
        for (a, A) in enumerate(AS)
            T = A == 0 ? T0 : zeros(o.nspec)
            if A != 0
                Threads.@threads for s in 1:o.nspec
                    T[s] = run_chain(S, base[s] .+ A .* sig[s] .* prof; compute_pvalues = false).T
                end
            end
            Ff[i, a, j] = mean(T .>= cutoff)
            Fm[i, a, j] = mean(mc_pvalue.(T, Ref(null)) .<= 0.05)
        end
        log_message("table (b): $(shape_name(sh)), w = $w: " * join([@sprintf("%.3f", Ff[i, a, j]) for a in eachindex(AS)], " "))
    end
    res = (; part = :b, shapes = shape_name.(SH), AS, WS, formula = Ff, mc = Fm, nspec = o.nspec, nnull = o.nnull,
           tau = o.tau, n = o.n, m = o.m, removed = n_removed_rfft(S.keep), kappa = S.kappa)
    save_result(result_path(o, "table_b" * shape_suffix(o)), res)
    write_table_csv(joinpath(o.outdir, "table_b" * shape_suffix(o) * ".csv"), res)
    return res
end

function write_table_csv(path, res)
    rows = Any[]
    for (i, sh) in enumerate(res.shapes), (a, A) in enumerate(res.AS), (j, w) in enumerate(res.WS)
        push!(rows, [sh, A, w, res.formula[i, a, j], res.mc[i, a, j]])
    end
    write_csv(path, ["shape", "A", "w", "f_formula", "f_mc"], permutedims(hcat(rows...)))
end

function merge_tables(paths::AbstractVector{<:AbstractString}, out::AbstractString)
    parts = load_result.(paths)
    shapes = unique(vcat([p.shapes for p in parts]...)); AS = parts[1].AS
    WS = sort(unique(vcat([p.WS for p in parts]...)))
    Ff = fill(NaN, length(shapes), length(AS), length(WS)); Fm = copy(Ff)
    for p in parts, (i, sh) in enumerate(p.shapes), (j, w) in enumerate(p.WS)
        I = findfirst(==(sh), shapes); Jw = findfirst(==(w), WS)
        Ff[I, :, Jw] = p.formula[i, :, j]; Fm[I, :, Jw] = p.mc[i, :, j]
    end
    res = merge(parts[1], (; shapes, WS, formula = Ff, mc = Fm))
    save_result(out, res)
    write_table_csv(replace(out, ".jls" => ".csv"), res)
    return res
end
