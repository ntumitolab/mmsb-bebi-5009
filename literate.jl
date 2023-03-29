using Distributed
using PrettyTables

@everywhere begin
    using Literate
end

folder = joinpath(@__DIR__, "docs")
nbs = [nb for nb in readdir(folder) if endswith(nb, ".jl")]

ts = pmap(nbs; on_error=ex->NaN) do nb
    @elapsed Literate.notebook(joinpath(folder, nb), folder; mdstrings=true)
end

pretty_table([nbs ts], header=["Notebook", "Elapsed (s)"])

any(isnan, ts) && error("Error(s) occured!")
