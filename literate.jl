using Distributed

# For all processes
@everywhere begin
    import Pkg
    Pkg.activate(@__DIR__)

    using Literate
    config = Dict("mdstrings" => true)
end

folder = joinpath(@__DIR__, "docs")
nbs = [nb for nb in readdir(folder) if splitext(nb)[end] == ".jl"]

ts = pmap(nbs; on_error=ex->NaN) do nb
    @elapsed Literate.notebook(joinpath(folder, nb), folder; config)
end

for (nb, t) in zip(nbs, ts)
    println(nb, " took ", t, " second(s). (0 means error occured)")
end

any(isnan, ts) && error("Error(s) occured!")
