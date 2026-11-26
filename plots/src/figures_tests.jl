function fig_nobkg_data(r)
    fig = Figure(size = (340, 300))
    a = Axis(fig[1, 1], xlabel = "bin", ylabel = L"x_i")
    lines!(a, 1:r.L, r.example; color = MIDGREY, linewidth = 0.4)
    i0, i1 = r.run
    w = max(1, i0 - 60):min(r.L, i1 + 60)
    vspan!(a, [i0], [i1]; color = (RED, 0.25))
    b = Axis(fig[2, 1], xlabel = "bin", ylabel = L"x_i", title = L"T_{\mathrm{obs}}=%$(round(r.T_example, digits = 2)),\ p=%$(round(r.p_example, digits = 3))")
    scatterlines!(b, collect(w), r.example[w]; color = MIDGREY, markersize = 3, linewidth = 0.6)
    scatterlines!(b, collect(i0:i1), r.example[i0:i1]; color = RED, markersize = 4, linewidth = 0.9)
    hlines!(b, [0]; color = :black, linewidth = 0.6)
    panel!(a, "a"); panel!(b, "b")
    return fig
end

function fig_nobkg_pvalues(r)
    fig = Figure(size = (340, 300))
    a = Axis(fig[1, 1], xlabel = L"p", ylabel = "density")
    phist!(a, p_hist(r.p_noise); color = :black)
    hlines!(a, [1]; color = MIDGREY, linestyle = :dash, linewidth = 0.8)
    b = Axis(fig[2, 1], xlabel = L"p", ylabel = "density", yscale = log10)
    phist!(b, p_hist(r.shapes["Gaussian"].p) .+ 1e-9; color = BLUE, label = L"A=20,\ w=4")
    phist!(b, p_hist(r.p_noise) .+ 1e-9; color = :black, label = "noise only")
    ylims!(b, 1e-2, 40)
    axislegend(b; position = :rt)
    panel!(a, "a"); panel!(b, "b")
    return fig
end
p_hist(p; nb = 20) = [count(x -> (k - 1) / nb <= x < k / nb || (k == nb && x == 1), p) for k in 1:nb]

function fig_shapes_pvalues(r)
    fig = Figure(size = (340, 300))
    names = ["Gaussian", "Lorentzian", "top-hat", "axion"]
    for (k, nm) in enumerate(names)
        haskey(r.shapes, nm) || continue
        ax = Axis(fig[div(k - 1, 2) + 1, (k - 1) % 2 + 1], xlabel = k > 2 ? L"p" : "", ylabel = isodd(k) ? "density" : "",
                  yscale = log10, title = nm)
        phist!(ax, p_hist(r.shapes[nm].p) .+ 1e-9; color = BLUE)
        phist!(ax, p_hist(r.p_noise) .+ 1e-9; color = :black)
        ylims!(ax, 1e-2, 40)
        ia = inset_axis(fig[div(k - 1, 2) + 1, (k - 1) % 2 + 1]; width = Relative(0.38), height = Relative(0.32))
        prof = r.shapes[nm].profile
        stairs!(ia, collect(eachindex(prof)) .- (length(prof) + 1) / 2, prof; step = :center, color = PURPLE)
        panel!(ax, string('a' + k - 1))
    end
    return fig
end

function fig_correlation(r)
    taus = r.figure_taus
    rows = [(lv, true) for lv in r.levels]; push!(rows, (r.level, false))
    nc = 1 + length(taus)
    fig = Figure(size = (40 * 2nc + 120, 40 * length(rows) + 110))
    cm = Reverse(:RdBu)
    for (side, tag) in enumerate(r.models), (i, (lv, bon)) in enumerate(rows)
        (bon || tag == "homo") || continue
        key = (tag, lv, bon)
        haskey(r.cases, key) || continue
        d = r.cases[key]
        for c in 1:nc
            M = c == 1 ? d.raw : d.taus[taus[c-1]].corr
            ax = Axis(fig[i, (side - 1) * nc + c], aspect = 1, yreversed = true,
                      title = i == 1 ? (c == 1 ? L"r" : tau_label(taus[c-1])) : "", titlesize = 8)
            hidedecorations!(ax); hidespines!(ax)
            heatmap!(ax, M; colormap = cm, colorrange = (-0.05, 0.05))
            c == 1 && Label(fig[i, (side - 1) * nc + c, Left()], bon ? level_label(lv) : L"b=0"; rotation = pi / 2,
                            fontsize = 8, padding = (0, 3, 0, 0))
        end
    end
    Label(fig[0, 1:nc], "homoscedastic"; fontsize = 10)
    length(r.models) > 1 && Label(fig[0, nc+1:2nc], "heteroscedastic"; fontsize = 10)
    Colorbar(fig[1:length(rows), 2nc + 1]; colormap = cm, limits = (-0.05, 0.05), label = "correlation")
    colgap!(fig.layout, 2); rowgap!(fig.layout, 2)
    return fig
end

function fig_corrdiff(r)
    fig = Figure(size = (720, 230))
    cm = Reverse(:RdBu); ncol = 0
    for (i, tag) in enumerate(r.models)
        t0 = tag == "homo" ? r.tau0.homo : r.tau0.hetero
        taus = sort(unique([r.figure_taus; r.tau_min; t0]))
        ncol = max(ncol, length(taus))
        d = r.cases[(tag, r.level, true)]
        for (c, tau) in enumerate(taus)
            e = d.taus[tau]
            ttl = tau == r.tau_min ? L"\tau_{\min}" : tau == t0 ? L"\tau_0" : tau_label(tau)
            ax = Axis(fig[i, c], aspect = 1, yreversed = true, title = i == 1 ? ttl : "", titlesize = 8)
            hidedecorations!(ax); hidespines!(ax)
            heatmap!(ax, e.diff; colormap = cm, colorrange = (-0.05, 0.05))
            text!(ax, 0.5, 0.06; text = @sprintf("%.4f", e.mean_diff), space = :relative, align = (:center, :bottom), fontsize = 7)
            c == 1 && Label(fig[i, c, Left()], tag; rotation = pi / 2, fontsize = 9, padding = (0, 3, 0, 0))
        end
    end
    Colorbar(fig[1:length(r.models), ncol + 1]; colormap = cm, limits = (-0.05, 0.05))
    colgap!(fig.layout, 3); rowgap!(fig.layout, 3)
    return fig
end

function fig_tauscan(r, t0 = nothing)
    fig = Figure(size = (560, 230))
    a = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = L"\tau", ylabel = L"\mathrm{rms}\,\rho_c(\ell),\ \ell=1\ldots250")
    b = Axis(fig[1, 2], xscale = log10, xlabel = L"\tau", ylabel = L"\mathrm{rms}\,\rho_c(\ell)")
    cols = palette(:viridis, length(r.levels))
    for (i, lv) in enumerate(r.levels), tag in r.models
        y = r.rms[(tag, lv)]
        ls = tag == "homo" ? :solid : :dash
        lines!(a, r.taus, pos(y); color = cols[i], linestyle = ls, label = tag == "homo" ? level_label(lv) : nothing)
        sel = (r.taus .>= r.tau_min / 20) .& (r.taus .<= r.tau_min * 20)
        scatter!(b, r.taus[sel], y[sel]; color = tag == "homo" ? cols[i] : :transparent, strokecolor = cols[i],
                 strokewidth = 0.8, markersize = 5)
    end
    lines!(a, r.taus, pos(r.rms_ref); color = :black, linestyle = :dot, label = "white noise")
    sel = (r.taus .>= r.tau_min / 20) .& (r.taus .<= r.tau_min * 20)
    lines!(b, r.taus[sel], r.rms_ref[sel]; color = :black, linestyle = :dot)
    for ax in (a, b)
        vlines!(ax, [r.tau_min]; color = RED, linestyle = :dash, linewidth = 0.9)
    end
    if t0 !== nothing
        for (tv, col) in ((t0.homo, GREEN), (t0.hetero, GREEN))
            vspan!(b, [tv / 1.04], [tv * 1.04]; color = (col, 0.3))
            vspan!(a, [tv / 1.04], [tv * 1.04]; color = (col, 0.3))
        end
    end
    axislegend(a; position = :rt, labelsize = 7.5)
    panel!(a, "a"); panel!(b, "b")
    return fig
end

function fig_fptp(r)
    fig = Figure(size = (720, 230))
    titles = [L"\mathcal{H}_0^{\mathrm{iid}}", L"\mathcal{H}_0^{\mathrm{filt}}", "boundary term removed"]
    sets = [[("iid_mc", ORANGE, "Monte Carlo"), ("pk", GREEN, L"P_{\mathcal{K}}"), ("pk_rm", VIOLET, L"P_{\mathcal{K}}\mathrm{,\ removed}")],
            [("filt_mc", BLUE, "Monte Carlo"), ("pk", GREEN, L"P_{\mathcal{K}}"), ("filt_mc_rm", LIGHTVIOLET, "MC, removed")],
            [("pk_rm", VIOLET, L"\mathcal{H}_0^{\mathrm{iid}}\ P_{\mathcal{K}}"), ("filt_mc_rm", LIGHTVIOLET, L"\mathcal{H}_0^{\mathrm{filt}}\ \mathrm{MC}")]]
    for (j, set) in enumerate(sets)
        ax = Axis(fig[1, j], xscale = log10, xlabel = L"\tau", ylabel = j == 1 ? L"f(p<0.05)" : "", title = titles[j])
        for (key, col, lab) in set, tag in keys(r.cases)
            d = r.cases[tag]; homo = tag == "homo"
            lines!(ax, r.taus, d.FP[key]; color = col, linestyle = homo ? :dash : :dot, linewidth = 1.1,
                   label = homo ? lab : nothing)
            scatter!(ax, r.taus, d.TP[key]; color = homo ? :transparent : col, strokecolor = col, strokewidth = 0.8,
                     markersize = 4.5)
        end
        hlines!(ax, [0.05]; color = MIDGREY, linestyle = :dot, linewidth = 0.8)
        vlines!(ax, [r.tau_min]; color = MIDGREY, linestyle = :dot, linewidth = 0.8)
        vlines!(ax, [0.05]; color = MIDGREY, linestyle = :dash, linewidth = 0.8)
        ylims!(ax, -0.02, 1.02)
        axislegend(ax; position = :lc, labelsize = 7.5)
        panel!(ax, string('a' + j - 1))
    end
    return fig
end

function fig_pvalues(r, null::Symbol, tag::AbstractString)
    fig = Figure(size = (300, 470))
    for (i, lv) in enumerate(r.levels)
        haskey(r.cases, (tag, lv)) || continue
        d = r.cases[(tag, lv)]
        pertau = hasproperty(d, :pertau) ? d.pertau : getproperty(d, Symbol("per\u03c4"))
        for (j, name) in enumerate(("noise", "signal"))
            ax = Axis(fig[i, j], yscale = log10, xlabel = i == length(r.levels) ? L"p" : "",
                      ylabel = j == 1 ? level_label(lv) : "", title = i == 1 ? (j == 1 ? "no signal" : "axion line") : "")
            ntau = length(d.taus)
            for (k, tau) in enumerate(d.taus)
                h = getfield(pertau[tau].hist[name], null)
                col = k == 1 ? RED : k == 2 ? GREEN : palette(:plasma, ntau - 2)[k-2]
                phist!(ax, h .+ 1e-9; color = col, linewidth = k <= 2 ? 1.3 : 0.7)
                if null == :iid
                    phist!(ax, getfield(pertau[tau].hist["b0"], null) .+ 1e-9; color = col, linestyle = :dash, linewidth = 0.6)
                end
            end
            phist!(ax, (name == "noise" ? d.raw_noise : d.raw_signal) .+ 1e-9; color = :black, linewidth = 1.0)
            hlines!(ax, [1]; color = MIDGREY, linestyle = :dash, linewidth = 0.7)
            ylims!(ax, 1e-2, 25)
            if j == 2
                ia = inset_axis(fig[i, j]; width = Relative(0.42), height = Relative(0.38))
                l = pertau[d.taus[end]].line
                lines!(ia, collect(eachindex(l)) .- div(length(l), 2), l; color = BLUE)
            end
            panel!(ax, string('a' + 2(i - 1) + j - 1))
        end
    end
    rowgap!(fig.layout, 3); colgap!(fig.layout, 4)
    return fig
end
