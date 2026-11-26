using MadmaxRuns, Printf
dirs = prepare_output_layout()
dir = isempty(ARGS) ? dirs.data : ARGS[1]
for part in ("a", "b")
    files = filter(f -> startswith(f, "table_$(part)_") && endswith(f, ".jls"), readdir(dir))
    isempty(files) && continue
    res = merge_tables(joinpath.(dir, files), joinpath(dir, "table_$part.jls"))
    println("table ($part): merged $(length(files)) parts")
end
include(joinpath(@__DIR__, "latex_table.jl"))
write_latex_table(dir)
