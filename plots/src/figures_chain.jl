function fig_examples(r)
    fig = Figure(size = (620, 860))
    tags = [t for t in ("homo", "hetero") if haskey(r.cases, t)]
    for (j, tag) in enumerate(tags)
        d = r.cases[tag]; nfull = length(d.x); m = div(nfull - length(d.r), 2)
        f = MadmaxRuns.PAPER_GRID.df_hz .* (0:nfull-1) ./ 1e6; fin = f[m+1:end-m]
        offset_label = L"f-f_{\mathrm{ref}}\ \mathrm{(MHz)}"
        a = Axis(fig[1, j], xlabel = offset_label, ylabel = j == 1 ? L"x" : "", title = tag == "homo" ? "homoscedastic" : "heteroscedastic")
        lines!(a, f, d.x; color = MIDGREY, linewidth = 0.4)
        lines!(a, fin, d.bhat; color = BLUE, linewidth = 1.0)
        c = argmax(d.line)
        scatter!(a, [f[c]], [maximum(d.x)]; marker = :dtriangle, color = RED, markersize = 8)
        ia = inset_axis(fig[1, j]; width = Relative(0.3), height = Relative(0.35), halign = 0.97, valign = 0.95)
        w = c-12:c+12
        scatter!(ia, f[w], d.x[w]; color = MIDGREY, markersize = 2.5)
        lines!(ia, f[w], d.b[w] .+ d.line[w]; color = RED)
        b = Axis(fig[2, j], xscale = log10, yscale = log10, xlabel = L"f\ \mathrm{(cycles/bin)}", ylabel = j == 1 ? "power" : "")
        k = 2:length(d.f)
        lines!(b, d.f[k], pos(d.Px[k]); color = LIGHTGREY, linewidth = 0.4, label = L"|\tilde{x}(f)|^2")
        lines!(b, d.f[k], pos(d.Pr[k]); color = BLUE, linewidth = 0.4, label = L"|\tilde{r}(f)|^2")
        lines!(b, d.f[k], pos(d.Pmodel[k]); color = :black, linestyle = :dash, label = L"\hat{\sigma}^2|G(f)|^2")
        hlines!(b, [mean(abs2, d.sigma_hat)]; color = RED, linestyle = :dot, label = L"\hat{\sigma}^2")
        j == 1 && axislegend(b; position = :rb, labelsize = 7)
        cc = Axis(fig[3, j], xlabel = offset_label, ylabel = j == 1 ? L"r_i,\ \bar\sigma r_i/\hat{\sigma}_i" : "")
        lines!(cc, fin, d.r; color = MIDGREY, linewidth = 0.3)
        lines!(cc, fin, d.pre; color = (BLUE, 0.6), linewidth = 0.3)
        text!(cc, 0.97, 0.05; space = :relative, align = (:right, :bottom), fontsize = 7,
              text = @sprintf("kurtosis %.2f \u2192 %.2f", d.kurt_r, d.kurt_pre))
        e = Axis(fig[4, j], xlabel = L"\langle b\rangle_w", ylabel = j == 1 ? L"\hat{\sigma}_w" : "")
        scatter!(e, d.bw, d.sw; color = (ORANGE, 0.6), markersize = 2.5)
        bb = range(extrema(d.bw)..., length = 40)
        lines!(e, bb, d.fit[1] .+ d.fit[2] .* bb; color = RED)
        lines!(e, bb, d.truth[1] .+ d.truth[2] .* bb; color = :black, linestyle = :dash)
        g = Axis(fig[5, j], xscale = log10, yscale = log10, xlabel = L"f\ \mathrm{(cycles/bin)}", ylabel = j == 1 ? "modes" : "")
        lines!(g, d.f[k], pos(d.before[k]); color = LIGHTGREY, linewidth = 0.4)
        lines!(g, d.f[k], pos(d.after[k]); color = BLUE, linewidth = 0.4)
        fmax = (maximum(d.removed) + 0.5) / (2 * (length(d.f) - 1))
        vspan!(g, [d.f[2]], [fmax]; color = (MIDGREY, 0.25))
        h = Axis(fig[6, j], yscale = log10, xlabel = L"p", ylabel = j == 1 ? "density" : "")
        phist!(h, p_hist(d.p_raw) .+ 1e-9; color = RED)
        phist!(h, p_hist(d.p_noise) .+ 1e-9; color = BLUE)
        phist!(h, p_hist(d.p_signal) .+ 1e-9; color = :black)
        hlines!(h, [1]; color = MIDGREY, linestyle = :dash, linewidth = 0.7)
        ylims!(h, 1e-2, 25)
        for (q, ax) in enumerate((a, b, cc, e, g, h))
            panel!(ax, string('a' + 6(j - 1) + q - 1))
        end
    end
    rowgap!(fig.layout, 4)
    return fig
end

function fig_leakage(r)
    fig = Figure(size = (340, 170))
    cols = palette(:viridis, length(r.ns))
    for (j, kind) in enumerate((:SG, :SGW))
        ax = Axis(fig[1, j], yscale = log10, xlabel = L"m", ylabel = j == 1 ? L"\phi" : "", title = string(kind))
        for (i, n) in enumerate(r.ns)
            scatterlines!(ax, r.ms, pos(r.phi[kind][i, :]); color = cols[i], markersize = 3, linewidth = 0.9, label = L"n=%$n")
        end
        hlines!(ax, [0.05]; color = MIDGREY, linestyle = :dot)
        panel!(ax, string('a' + j - 1))
        j == 2 && axislegend(ax; position = :rb, labelsize = 7)
    end
    return fig
end

function fig_tolerance(r)
    fig = Figure(size = (620, 330))
    f = r.f ./ 1e6
    a = Axis(fig[1, 1], yscale = log10, xlabel = L"f\ \mathrm{(MHz)}", ylabel = L"\Delta_i")
    lines!(a, f, pos(r.Delta); color = BLUE, label = "exact")
    lines!(a, f, pos(r.Delta_taylor); color = RED, linestyle = :dash, label = "Eq. (t4)")
    axislegend(a; position = :rt)
    b = Axis(fig[1, 2], yscale = log10, xlabel = L"f\ \mathrm{(MHz)}", ylabel = L"|\varepsilon_i|")
    lines!(b, f, pos(r.eps_max); color = BLUE, label = "Eq. (t2)")
    lines!(b, f, pos(r.eps_taylor_max); color = RED, linestyle = :dash, label = "Eq. (t5)")
    lines!(b, f, pos(r.bound_max); color = :black, linestyle = :dot, label = "bound (t3)")
    axislegend(b; position = :rt)
    cols = palette(:viridis, length(r.levels))
    c = Axis(fig[2, 1], yscale = log10, xlabel = L"f\ \mathrm{(MHz)}", ylabel = L"|(G\mathbf{b})_i|/\sigma_i")
    d = Axis(fig[2, 2], yscale = log10, xlabel = L"f\ \mathrm{(MHz)}", ylabel = "sliding mean ratio")
    for (i, lv) in enumerate(r.levels)
        lines!(c, f, pos(r.Gb_sigma_max[lv]); color = cols[i], label = level_label(lv))
        lines!(d, f, pos(r.sliding_max[lv]); color = cols[i])
    end
    axislegend(c; position = :rt, labelsize = 7)
    for (q, ax) in enumerate((a, b, c, d))
        panel!(ax, string('a' + q - 1))
    end
    return fig
end

function fig_tau0_map(r)
    fig = Figure(size = (300, 230))
    ax = Axis(fig[1, 1], yscale = log10, xlabel = L"m", ylabel = L"\tau_0")
    cols = palette(:viridis, length(r.ns))
    for (i, n) in enumerate(r.ns)
        scatterlines!(ax, r.ms, pos(r.map["homo"][i, :]); color = cols[i], markercolor = :transparent, strokecolor = cols[i],
                      strokewidth = 0.8, markersize = 4, linewidth = 0.8, label = L"n=%$n")
        scatterlines!(ax, r.ms, pos(r.map["hetero"][i, :]); color = cols[i], markersize = 4, linewidth = 0.8, linestyle = :dash)
    end
    if haskey(r.tau0, (r.n, r.m, "homo"))
        scatter!(ax, [r.m, r.m], [r.tau0[(r.n, r.m, "homo")], r.tau0[(r.n, r.m, "hetero")]]; marker = :rect,
                 color = :transparent, strokecolor = :black, strokewidth = 1.2, markersize = 9)
    end
    Legend(fig[1, 2], ax; labelsize = 7.5)
    return fig
end

function fig_tau0_crossing(r)
    fig = Figure(size = (300, 230))
    ax = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = L"\tau", ylabel = L"\beta_{\max}\sqrt{\ell_\tau}")
    for (tag, col) in (("homo", :black), ("hetero", RED))
        d = r.cross[tag]
        lines!(ax, r.taus, pos(d.significance); color = col, label = tag)
        lines!(ax, r.taus, pos(d.closed); color = col, linestyle = :dash)
        vlines!(ax, [d.tau0]; color = col, linestyle = :dot, linewidth = 0.8)
    end
    lines!(ax, r.taus, r.cross["homo"].z_alpha; color = GREEN, label = L"z_\alpha")
    vlines!(ax, [r.tau_min]; color = RED, linestyle = :dash, linewidth = 0.8)
    axislegend(ax; position = :rt)
    return fig
end

function fig_boundary(r)
    fig = Figure(size = (520, 300))
    a = Axis(fig[1, 1], ylabel = L"z_i")
    lines!(a, r.rows, r.z; color = MIDGREY, linewidth = 0.3)
    lines!(a, r.rows, r.component; color = RED, linewidth = 1.2)
    b = Axis(fig[2, 1], xlabel = "bin", ylabel = L"z'_i")
    lines!(b, r.rows, r.zp; color = BLUE, linewidth = 0.3)
    text!(a, 0.97, 0.95; space = :relative, align = (:right, :top), fontsize = 8, text = @sprintf("T = %.1f", r.T_before))
    text!(b, 0.97, 0.95; space = :relative, align = (:right, :top), fontsize = 8, text = @sprintf("T = %.1f", r.T_after))
    linkxaxes!(a, b); hidexdecorations!(a; ticks = false, minorticks = false)
    panel!(a, "a"); panel!(b, "b")
    return fig
end

function fig_boundary_svd(r)
    fig = Figure(size = (340, 170))
    a = Axis(fig[1, 1], yscale = log10, xlabel = L"k", ylabel = L"s_k")
    k = 1:min(20, length(r.s))
    scatterlines!(a, k, pos(r.s[k]); color = BLUE, markersize = 4)
    b = Axis(fig[1, 2], yscale = log10, xlabel = L"k", ylabel = "uncaptured variance")
    scatterlines!(b, k, pos(r.uncaptured[k]); color = RED, markersize = 4)
    vlines!(b, [6]; color = MIDGREY, linestyle = :dot)
    panel!(a, "a"); panel!(b, "b")
    return fig
end

function fig_tcut(r)
    fig = Figure(size = (300, 400))
    a = Axis(fig[1, 1], xlabel = L"T", ylabel = "density")
    b = Axis(fig[2, 1], xlabel = L"T", ylabel = L"P(T'\geq T)", yscale = log10)
    cols = palette(:viridis, length(r.ntau))
    for (i, nt) in enumerate(r.ntau)
        d = r.cases[nt]
        x = (d.edges[1:end-1] .+ d.edges[2:end]) ./ 2
        lines!(a, x, d.density; color = cols[i], label = L"n_\tau=%$nt")
        lines!(a, mids(r.tgrid), density_from_survival(r.tgrid, d.S_exact); color = cols[i], linestyle = :dash)
        lines!(a, mids(r.tgrid), density_from_survival(r.tgrid, d.S_first); color = cols[i], linestyle = :dot)
        lines!(b, r.tgrid, pos(d.S_mc); color = cols[i])
        lines!(b, r.tgrid, pos(d.S_exact); color = cols[i], linestyle = :dash)
        lines!(b, r.tgrid, pos(d.S_first); color = cols[i], linestyle = :dot)
    end
    d1 = r.cases[r.ntau[1]]
    lines!(a, mids(r.tgrid), density_from_survival(r.tgrid, d1.S_pbc); color = :black)
    lines!(b, r.tgrid, pos(d1.S_pbc); color = :black, label = "no modes removed")
    axislegend(a; position = :rt, labelsize = 7.5)
    panel!(a, "a"); panel!(b, "b")
    return fig
end
mids(t) = (t[1:end-1] .+ t[2:end]) ./ 2
density_from_survival(t, S) = -diff(S) ./ diff(t)

function fig_g2nd(r)
    fig = Figure(size = (300, 420))
    for (i, nt) in enumerate(r.ntau)
        d = r.cases[nt]
        ax = Axis(fig[i, 1], xlabel = i == length(r.ntau) ? L"t" : "", ylabel = L"P(T\geq t)/P_{\mathrm{MC}}", title = L"n_\tau=%$nt")
        ok = d.S_mc .> 0
        band!(ax, r.tgrid[ok], 1 .- d.band[ok], 1 .+ d.band[ok]; color = (LIGHTGREY, 0.6))
        lines!(ax, r.tgrid[ok], d.ratio_exact[ok]; color = BLUE, label = "exact")
        lines!(ax, r.tgrid[ok], d.ratio_second[ok]; color = GREEN, label = "second order")
        lines!(ax, r.tgrid[ok], d.ratio_first[ok]; color = RED, linestyle = :dash, label = "first order")
        hlines!(ax, [1]; color = :black, linewidth = 0.6)
        i == 1 && axislegend(ax; position = :lb, labelsize = 7.5)
        panel!(ax, string('a' + i - 1))
    end
    return fig
end
