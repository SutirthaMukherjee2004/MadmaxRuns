function example_panel(layout, position; kw...)
    return Axis(layout[position...]; xlabelsize = 6.5, ylabelsize = 6.5, titlesize = 6.5,
                xticklabelsize = 5.5, yticklabelsize = 5.5, xticksize = 2.5, yticksize = 2.5,
                xminorticksize = 1.2, yminorticksize = 1.2, kw...)
end

function example_residual_density!(ax, axis, values, color)
    lo, hi = axis
    nb = 50
    edges = range(lo, hi; length = nb + 1)
    h = zeros(Int, nb)
    for x in values
        lo <= x <= hi || continue
        k = clamp(searchsortedlast(edges, x), 1, nb)
        h[k] += 1
    end
    counts = repeat(h; inner = 2)
    ys = collect(Iterators.flatten((edges[k], edges[k+1]) for k in 1:nb))
    lines!(ax, counts, ys; color, linewidth = 0.6)
end

function fig_examples(r)
    tags = [tag for tag in ("homo", "hetero") if haskey(r.cases, tag)]
    fig = Figure(size = (length(tags) == 2 ? 456 : 228, 428), figure_padding = 6)
    if length(tags) == 2
        Box(fig[1, 2]; color = :black, strokewidth = 0, width = 0.55, tellwidth = false)
        colsize!(fig.layout, 2, Fixed(8))
    end
    for (j, tag) in enumerate(tags)
        grid = GridLayout(fig[1, 2j-1])
        d = r.cases[tag]
        nfull = length(d.x)
        m = div(nfull - length(d.r), 2)
        f = MadmaxRuns.PAPER_GRID.df_hz .* (0:nfull-1) ./ 1e6
        fin = f[m+1:end-m]
        model = tag == "homo" ? "homoscedastic" : "heteroscedastic"
        level = @sprintf("%g", 100r.level)
        offset_label = L"f-f_{\mathrm{ref}}\ \mathrm{[MHz]}"
        mode_label = L"\mathrm{mode}\ f\ \mathrm{[1/bin]}"
        title = L"%$model,\ \sigma/b_{\mathrm{med}}=%$level\%"

        a = example_panel(grid, (1, 1:2); xlabel = offset_label, ylabel = L"x\ [10^{-22}\,\mathrm{W}]",
                          title, yticks = 0:25:175)
        lines!(a, f, d.x; color = MIDGREY, linewidth = 0.28, label = L"x=b+n+s")
        lines!(a, fin, d.bhat; color = BLUE, linewidth = 0.65, label = L"\hat{b}=Sx")
        xlims!(a, -1, 51)
        ylims!(a, 0, 175)
        axislegend(a; position = :rt, labelsize = 4.8, patchsize = (9, 5))
        c = argmax(d.line)
        scatter!(a, [f[c]], [100]; marker = :dtriangle, color = RED, markersize = 6)
        text!(a, f[c], 127; text = L"\mathrm{axion\ line\ injected}\ (3\sigma_i)", color = RED,
              align = (:center, :bottom), fontsize = 5)
        lines!(a, [f[c], f[c]], [102, 120]; color = RED, linewidth = 0.65)
        lo, hi = f[c] - 0.25, f[c] + 0.5
        w = searchsortedfirst(f, lo):searchsortedlast(f, hi)
        ia = inset_axis(grid[1, 1:2]; width = Relative(0.33), height = Relative(0.36),
                        halign = 0.92, valign = 0.25, xticks = [10.25, 10.5, 10.75, 11.0],
                        xticklabelsize = 4.2, yticklabelsize = 4.2)
        lines!(ia, f[w], d.x[w]; color = GREY, linewidth = 0.24)
        lines!(ia, f[w], d.b[w] .+ d.line[w]; color = RED, linewidth = 0.6)
        xlims!(ia, lo, hi)

        xp = example_panel(grid, (2, 1); xlabel = mode_label, xscale = log10, yscale = log10,
                           title = L"|\tilde{x}(f)|^2:\ \mathrm{background\ at\ low}\ f")
        rp = example_panel(grid, (2, 2); xlabel = mode_label, xscale = log10, yscale = log10,
                           title = L"|\tilde{r}(f)|^2\ \simeq\ \hat{\sigma}^2 |G|^2")
        k = 2:length(d.f)
        for ax in (xp, rp)
            xlims!(ax, 1e-4, 0.55)
        end
        lines!(xp, d.f[k], pos(d.Px[k]); color = BLUE, linewidth = 0.28)
        hlines!(xp, [mean(abs2, d.sigma_hat)]; color = RED, linestyle = :dash, linewidth = 0.55)
        ylims!(xp, 1e-3, 2e7)
        lines!(rp, d.f[k], pos(d.Pr[k]); color = BLUE, linewidth = 0.28)
        lines!(rp, d.f[k], pos(d.Pmodel[k]); color = RED, linestyle = :dash, linewidth = 0.55)
        ylims!(rp, 1e-5, 1e3)

        residual = GridLayout(grid[3, 1])
        cc = example_panel(residual, (1, 1); xlabel = offset_label,
                           title = L"r\ \mathrm{(grey)},\ \bar{\sigma}r_i/\hat{\sigma}_i\ \mathrm{(blue)}")
        dd = example_panel(residual, (1, 2))
        colsize!(residual, 2, Relative(0.18))
        colgap!(residual, 0)
        hidedecorations!(dd)
        lines!(cc, fin, d.r; color = (MIDGREY, 0.62), linewidth = 0.22)
        lines!(cc, fin, d.pre; color = (BLUE, 0.5), linewidth = 0.22)
        lines!(cc, fin, 2 .* d.sigma_hat; color = RED, linewidth = 0.5)
        lines!(cc, fin, -2 .* d.sigma_hat; color = RED, linewidth = 0.5)
        limit = max(18.0, 5ceil(maximum(abs, d.r) / 5))
        rlimits = (-limit, limit)
        xlims!(cc, 0, 50)
        ylims!(cc, rlimits...)
        ylims!(dd, rlimits...)
        example_residual_density!(dd, rlimits, d.r, MIDGREY)
        example_residual_density!(dd, rlimits, d.pre, BLUE)

        e = example_panel(grid, (3, 2); xlabel = L"\langle b\rangle_w",
                          title = L"\mathrm{fit}\ \hat{\sigma}_w=\hat{\sigma}_0+\hat{\sigma}_1\langle b\rangle_w")
        scatter!(e, d.bw, d.sw; color = ORANGE, markersize = 1.7, label = L"\mathrm{windows},\ W=96")
        bb = range(0, max(105, maximum(d.bw)); length = 40)
        lines!(e, bb, d.fit[1] .+ d.fit[2] .* bb; color = RED, linewidth = 0.65, label = "fit")
        lines!(e, bb, d.truth[1] .+ d.truth[2] .* bb; color = :black, linestyle = :dash, linewidth = 0.65, label = "truth")
        xlims!(e, first(bb), last(bb))
        ylims!(e, 0, max(4, 1.04maximum(d.sw)))
        axislegend(e; position = :rb, labelsize = 4.8, patchsize = (9, 5))

        g = example_panel(grid, (4, 1); xlabel = mode_label, xscale = log10, yscale = log10,
                         title = L"\mathrm{before/after}\ \div G:\ |\tilde{c}|^2\ \mathrm{on}\ M")
        removed = BitSet(d.removed)
        first_retained = findfirst(i -> !in(i - 1, removed), eachindex(d.f))
        if first_retained !== nothing && d.f[first_retained] > 1e-4
            vspan!(g, [1e-4], [d.f[first_retained]]; color = (MIDGREY, 0.28))
        end
        lines!(g, d.f[k], pos(d.before[k]); color = LIGHTGREY, linewidth = 0.28)
        lines!(g, d.f[k], pos(d.after[k]); color = (BLUE, 0.78), linewidth = 0.28)
        hlines!(g, [1.0]; color = RED, linestyle = :dash, linewidth = 0.55)
        xlims!(g, 1e-4, 0.55)
        ylims!(g, 1e-4, 350)

        h = example_panel(grid, (4, 2); yscale = log10, xlabel = L"p", ylabel = "counts",
                          title = L"p\ \mathrm{values},\ H_0^{\mathrm{iid}}\ \mathrm{(bin\ order)}")
        raw = p_hist(d.p_raw)
        edges = range(0, 1; length = length(raw) + 1)
        mid = (edges[1:end-1] .+ edges[2:end]) ./ 2
        barplot!(h, mid, raw .+ 1e-9; width = step(edges), fillto = 1, gap = 0,
                 color = (RED, 0.35), strokewidth = 0, label = L"\mathrm{raw}\ r,\ \mathrm{no\ signal}")
        phist!(h, p_hist(d.p_noise) .+ 1e-9; counts = true, color = BLUE, label = "whitened, no signal")
        phist!(h, p_hist(d.p_signal) .+ 1e-9; counts = true, color = :black, label = "whitened, axion injected")
        xlims!(h, 0, 1)
        ylims!(h, 0.8, 1e4)
        axislegend(h; position = :lt, labelsize = 4.2, patchsize = (9, 5))
        colgap!(grid, 17)
        rowgap!(grid, 14)
    end
    colgap!(fig.layout, 9)
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
