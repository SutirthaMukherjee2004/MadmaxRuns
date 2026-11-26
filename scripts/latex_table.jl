using MadmaxRuns, Printf

function write_latex_table(dir::AbstractString)
    parts = [(p, joinpath(dir, "table_$p.jls")) for p in ("a", "b")]
    parts = filter(x -> isfile(x[2]), parts)
    isempty(parts) && return println("no table results in $dir")
    io = IOBuffer()
    first = load_result(parts[1][2])
    AS = first.AS
    println(io, "\\begin{tabular}{ll", "c"^length(AS), "}")
    println(io, "\\toprule")
    println(io, "shape & \$w\$ & ", join(["\$A=$(Int(A))\$" for A in AS], " & "), " \\\\")
    for (p, path) in parts
        r = load_result(path)
        println(io, "\\midrule")
        println(io, "\\multicolumn{$(2 + length(AS))}{l}{(", p, ") ", p == "a" ? "no background, \$P_{\\rm BC}\$" :
                    "background, \$P_{\\mathcal K}\$", "} \\\\")
        for (i, sh) in enumerate(r.shapes), (j, w) in enumerate(r.WS)
            cells = [isnan(r.formula[i, a, j]) ? "--" : @sprintf("%.3f (%.3f)", r.formula[i, a, j], r.mc[i, a, j]) for a in eachindex(r.AS)]
            println(io, j == 1 ? sh : "", " & ", @sprintf("%g", w), " & ", join(cells, " & "), " \\\\")
        end
    end
    println(io, "\\bottomrule\n\\end{tabular}")
    path = joinpath(dir, "table_shapes4.tex")
    write(path, take!(io))
    println("wrote $path")
end
