const BLUE = Makie.to_color("#1b6ca8")
const RED = Makie.to_color("#b3123b")
const ORANGE = Makie.to_color("#c2410c")
const GREEN = Makie.to_color("#2e8b57")
const PURPLE = Makie.to_color("#7030A0")
const VIOLET = Makie.to_color("#6a3d9a")
const LIGHTVIOLET = Makie.to_color("#b39ddb")
const GREY = Makie.to_color("#606060")
const MIDGREY = Makie.to_color("#8c8c8c")
const LIGHTGREY = Makie.to_color("#b5b5b5")

function serif_fonts()
    fonts = (; regular = "DejaVu Serif", bold = "DejaVu Serif Bold", italic = "DejaVu Serif Italic")
    FT = Makie.FreeTypeAbstraction
    for (style, name) in pairs(fonts)
        face = FT.findfont(name)
        face === nothing && return nothing
        FT.family_name(face) == "dejavu serif" || return nothing
        if style != :regular
            occursin(string(style), FT.style_name(face)) || return nothing
        end
    end
    return fonts
end

function paper_theme()
    tick = (; xtickalign = 1, ytickalign = 1, xticksmirrored = true, yticksmirrored = true,
            xminorticksvisible = true, yminorticksvisible = true, xminortickalign = 1, yminortickalign = 1,
            xticksize = 3.5, yticksize = 3.5, xminorticksize = 2, yminorticksize = 2,
            xtickwidth = 0.7, ytickwidth = 0.7, xminortickwidth = 0.6, yminortickwidth = 0.6,
            xgridvisible = false, ygridvisible = false, spinewidth = 0.8,
            xlabelsize = 10, ylabelsize = 10, xticklabelsize = 9, yticklabelsize = 9, titlesize = 10,
            titlefont = :regular)
    fonts = serif_fonts()
    th = Theme(fontsize = 9.5, figure_padding = 4,
                 Axis = tick,
                 Legend = (; framevisible = false, labelsize = 8.5, patchsize = (14, 8), rowgap = 0, padding = (2, 2, 2, 2)),
                 Colorbar = (; ticklabelsize = 8.5, labelsize = 9.5, tickalign = 1, ticksize = 3, spinewidth = 0.7),
                 Lines = (; linewidth = 1.2), Scatter = (; markersize = 5), Stairs = (; linewidth = 1.1))
    return fonts === nothing ? th : merge(Theme(fonts = fonts), th)
end

panel!(ax, s::AbstractString) = text!(ax, 0.035, 0.965; text = "($s)", space = :relative, align = (:left, :top),
                                      fontsize = 10, font = :regular)

palette(cmap::Symbol, n::Integer) = n == 1 ? [cgrad(cmap)[0.0]] : [cgrad(cmap)[x] for x in range(0, 0.85, length = n)]

function inset_axis(pos; width = Relative(0.36), height = Relative(0.36), halign = 0.95, valign = 0.92, kw...)
    ax = Axis(pos; width, height, halign, valign, xticklabelsize = 7, yticklabelsize = 7,
              xticksize = 2, yticksize = 2, xminorticksvisible = false, yminorticksvisible = false,
              backgroundcolor = :white, kw...)
    translate!(ax.blockscene, 0, 0, 150)
    return ax
end

pos(v) = [x > 0 && isfinite(x) ? Float64(x) : NaN for x in v]

level_label(lv) = lv < 0.01 ? @sprintf("%.2f %%", 100lv) : @sprintf("%g %%", round(100lv, sigdigits = 3))
tau_label(tau) = (e = floor(Int, log10(tau)); m = tau / 10.0^e; isapprox(m, 1; atol = 0.05) ? L"\tau=10^{%$e}" :
                L"\tau=%$(round(m, digits = 1))\times10^{%$e}")

function phist!(ax, h::AbstractVector; counts::Bool = false, kw...)
    nb = length(h); e = range(0, 1, length = nb + 1)
    d = counts ? h : h ./ (sum(h) / nb)
    stairs!(ax, collect(e), [d; d[end]]; step = :post, kw...)
end

function savefig(fig, dir, name)
    mkpath(dir)
    pdf = name * ".pdf"
    png = name * ".png"
    mktempdir(dir) do staging
        save(joinpath(staging, pdf), fig; pt_per_unit = 1)
        save(joinpath(staging, png), fig; px_per_unit = 300 / 72)
        for file in (pdf, png)
            Base.Filesystem.rename(joinpath(staging, file), joinpath(dir, file))
        end
    end
    return joinpath(dir, pdf)
end
