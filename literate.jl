using PrettyTables
using Literate

basedir = "docs/"
config = Dict("mdstrings" => true, "execute" => true)

nbs = String[]
ts = Float64[]

# Collect the list of Literate notebooks
for (root, dirs, files) in walkdir(basedir)
    for file in files
        if (endswith(file, ".jl"))
            nb = joinpath(root, file)
            push!(nbs, nb)
        end
    end
end

# Execute notebooks
for nb in nbs
    try
        t = withenv("JULIA_DEBUG" => "Literate") do
            @elapsed Literate.notebook(nb, dirname(nb); config)
        end
        push!(ts, t)
    catch e
        println("An error occured in the cell above:", e)
        push!(ts, NaN)
    end
end

pretty_table([nbs ts], header=["Notebook", "Elapsed (s)"])

any(isnan, ts) && error("Please checkout errors.")
