function fig_filter(r)
    fig = Figure(size = (340, 250))
    a = Axis(fig[1, 1], xlabel = L"d\ \mathrm{(bins)}", ylabel = L"S_d")
    b = Axis(fig[1, 2], xlabel = L"f\ \mathrm{(cycles/bin)}", ylabel = L"|G(f)|^2", xscale = log10, yscale = log10)
    cols = [BLUE, RED, GREEN, ORANGE]
    for (k, F) in enumerate(r.filters)
        lab = "$(F.kind) ($(F.n), $(F.m))"
        lines!(a, F.d, F.S; color = cols[k], label = lab)
        keep = 2:length(F.f)
        lines!(b, F.f[keep], pos(F.G2[keep]); color = cols[k])
    end
    ylims!(b, 1e-14, 3)
    panel!(a, "a"); panel!(b, "b")
    axislegend(a; position = :rt)
    return fig
end

function fig_iteration(r)
    fig = Figure(size = (340, 170))
    lab = if hasproperty(r, :relative_changes) && r.relative_changes
        [L"\mathrm{median}\ |\hat{\sigma}_1^{(k)}/\hat{\sigma}_1^{(k-1)}-1|",
         L"\mathrm{median}\ |\hat{\sigma}_0^{(k)}/\hat{\sigma}_0^{(k-1)}-1|"]
    else
        [L"\mathrm{median}\ |\Delta\hat{\sigma}_1^{(k)}|", L"\mathrm{median}\ |\Delta\hat{\sigma}_0^{(k)}|"]
    end
    axs = [Axis(fig[1, j], xlabel = L"k", ylabel = lab[j], yscale = log10) for j in 1:2]
    cols = palette(:viridis, length(r.levels))
    for (i, lv) in enumerate(r.levels), (j, key) in enumerate((:dsigma1, :dsigma0))
        k = 1:r.nsteps
        l = getfield(r.cases[(:lsq, lv)], key); w = getfield(r.cases[(:window, lv)], key)
        scatterlines!(axs[j], k, pos(l); color = cols[i], marker = :circle, markersize = 5, linewidth = 0.9,
                      label = level_label(lv))
        scatterlines!(axs[j], k, pos(w); color = cols[i], marker = :rect, markercolor = :transparent,
                      strokecolor = cols[i], strokewidth = 0.9, markersize = 5, linewidth = 0.9, linestyle = :dash)
    end
    panel!(axs[1], "a"); panel!(axs[2], "b")
    axislegend(axs[1]; position = :rt)
    return fig
end

function fig_background(r)
    fig = Figure(size = (500, 420))
    reference = first(r.f)
    MHz = (r.f .- reference) ./ 1e6
    offset_label = L"f-f_{\mathrm{ref}}\ \mathrm{(MHz)}"
    a = Axis(fig[1, 1:2], xlabel = offset_label, ylabel = L"b\ (10^{-22}\,\mathrm{W})")
    for b in r.reals
        lines!(a, MHz, b; color = LIGHTGREY, linewidth = 0.6)
    end
    lines!(a, MHz, r.nominal; color = BLUE, linewidth = 1.2)
    ia = inset_axis(fig[1, 1:2]; width = Relative(0.25), height = Relative(0.35), halign = 0.97, valign = 0.95)
    c = round(Int, r.centre)
    w = c-15:c+15
    stairs!(ia, MHz[w], r.lines[1][w] ./ r.sigma[w]; step = :center, color = PURPLE)
    hidedecorations!(ia; ticklabels = false, ticks = false)
    b = Axis(fig[2, 1], xlabel = offset_label, ylabel = L"x\ (10^{-22}\,\mathrm{W})")
    z = c-250:c+250
    lines!(b, MHz[z], r.x[z]; color = MIDGREY, linewidth = 0.6, label = L"b+n")
    lines!(b, MHz[z], r.nominal[z]; color = BLUE, label = L"b")
    for (k, (l, col)) in enumerate(zip(r.lines, palette(:viridis, 3)))
        lines!(b, MHz[z], r.nominal[z] .+ l[z]; color = col, label = L"%$(Int(round(r.widths[k] / r.widths[1])))\times w_{\mathrm{SHM}}")
    end
    axislegend(b; position = :lt, labelsize = 7.5)
    cax = Axis(fig[2, 2], xlabel = offset_label, ylabel = L"r=G\mathbf{x}")
    fr = (r.f_r .- reference) ./ 1e6
    lines!(cax, fr, r.r; color = MIDGREY, linewidth = 0.4)
    lines!(cax, fr, r.envelope; color = RED, linewidth = 1.0)
    lines!(cax, fr, -r.envelope; color = RED, linewidth = 1.0)
    panel!(a, "a"); panel!(b, "b"); panel!(cax, "c")
    return fig
end

function fig_spectra(r)
    fig = Figure(size = (720, 300))
    for (i, tag) in enumerate(r.models), (j, lv) in enumerate(r.levels)
        d = r.spectra[(tag, lv)]
        ax = Axis(fig[i, j], xlabel = i == length(r.models) ? L"f\ \mathrm{(MHz)}" : "",
                  ylabel = j == 1 ? L"x\ (10^{-22}\,\mathrm{W})" : "",
                  title = "$(tag == "homo" ? "homoscedastic" : "heteroscedastic"), " * level_label(lv))
        lines!(ax, r.f ./ 1e6, d.x; color = MIDGREY, linewidth = 0.4)
        lines!(ax, r.f ./ 1e6, d.b; color = BLUE, linewidth = 1.0)
        panel!(ax, string('a' + (i - 1) * length(r.levels) + j - 1))
    end
    return fig
end

function fig_sigmafit(r)
    fig = Figure(size = (720, 330))
    for (i, tag) in enumerate(r.models), (j, lv) in enumerate(r.levels)
        d = r.cases[(tag, lv)]
        ax = Axis(fig[i, j], xlabel = i == length(r.models) ? L"\langle b\rangle_w" : "",
                  ylabel = j == 1 ? L"\hat{\sigma}_w" : "",
                  title = "$(tag), " * level_label(lv))
        scatter!(ax, d.bw256, d.sw256; color = (BLUE, 0.6), markersize = 3, label = L"W=256")
        scatter!(ax, d.bw96, d.sw96; color = (ORANGE, 0.5), markersize = 2.5, label = L"W=96")
        scatter!(ax, d.bblk, d.brms; color = (PURPLE, 0.6), marker = :utriangle, markersize = 3, label = "block rms")
        bb = range(extrema(d.bw96)..., length = 50)
        lines!(ax, bb, d.window[1] .+ d.window[2] .* bb; color = RED, label = "window fit")
        lines!(ax, bb, d.lsq[1] .+ d.lsq[2] .* bb; color = PURPLE, label = "Eq. (s6)")
        lines!(ax, bb, d.truth[1] .+ d.truth[2] .* bb; color = :black, linestyle = :dash, label = "truth")
        vlines!(ax, [r.bmed]; color = MIDGREY, linestyle = :dot, linewidth = 0.8)
        panel!(ax, string('a' + (i - 1) * length(r.levels) + j - 1))
        (i, j) == (1, 1) && axislegend(ax; position = :lt, labelsize = 7)
    end
    return fig
end
