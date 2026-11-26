const NOBKG_L = 10_000

function study_nobkg(o::StudyOptions)
    L = NOBKG_L; cen = L / 2
    ex = randn(spectrum_rng(o, 0, 1), L)
    Tex = squares_T(ex); _, i0, i1 = squares_T_run(ex)
    tab = pbc_table(L; tmax = 120.0)
    log_message(@sprintf("nobkg: example T_obs = %.2f in bins %d\u2013%d, p = %.4f", Tex, i0, i1, pbc_pvalue(Tex, L)))
    T0 = zeros(o.nspec)
    Threads.@threads for s in 1:o.nspec
        T0[s] = squares_T(randn(spectrum_rng(o, 1, s), L))
    end
    p0 = pvalue.(Ref(tab), T0)
    shapes = Dict{String,Any}()
    for sh in ALL_SHAPES
        prof = line_profile(sh, 20.0, 4.0, L, cen)
        T = zeros(o.nspec)
        Threads.@threads for s in 1:o.nspec
            T[s] = squares_T(randn(spectrum_rng(o, 1, s), L) .+ prof)
        end
        p = pvalue.(Ref(tab), T)
        shapes[shape_name(sh)] = (; p, f05 = fraction_below(p, 0.05), profile = prof[Int(cen)-30:Int(cen)+30])
        log_message(@sprintf("nobkg: %-10s A = 20, w = 4: f(p<=0.05) = %.3f", shape_name(sh), fraction_below(p, 0.05)))
    end
    Tn = zeros(o.nnull)
    Threads.@threads for s in 1:o.nnull
        Tn[s] = squares_T(randn(spectrum_rng(o, 2, s), L))
    end
    null = NullSet(:iid, sort!(Tn))
    pmc = mc_pvalue.(T0, Ref(null))
    log_message(@sprintf("nobkg: noise only f(p<=0.05) = %.3f (P_BC), %.3f (MC)", fraction_below(p0, 0.05), fraction_below(pmc, 0.05)))
    res = (; L, example = ex, T_example = Tex, run = (i0, i1), p_example = pbc_pvalue(Tex, L),
           p_noise = p0, p_noise_mc = pmc, shapes)
    save_result(result_path(o, "nobkg"), res)
    return res
end
