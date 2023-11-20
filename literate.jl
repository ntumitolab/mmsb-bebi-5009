using Distributed
using PrettyTables

@everywhere begin
    using Literate
end

basedir = "docs"
config = Dict("mdstrings" => true, "execute" => true, "flavor" => Literate.CommonMarkFlavor())

nbs = String[]

# Collect the list of Literate notebooks (ends with .jl)
for (root, dirs, files) in walkdir(basedir)
    for file in files
        if (endswith(file, ".jl"))
            push!(nbs, joinpath(root, file))
        end
    end
end

# Execute the notebooks in worker processes
ts = pmap(nbs; on_error=ex->NaN) do nb
    cd(dirname(nb)) do
        @elapsed Literate.markdown(basename(nb); config)
    end
end

pretty_table([nbs ts], header=["Notebook", "Elapsed (s)"])

# Debug notebooks one by one if there are errors
for (nb, t) in zip(nbs, ts)
    if isnan(t)
        println("Debugging notebook: ", nb)
        try
            withenv("JULIA_DEBUG" => "Literate") do
                cd(dirname(nb)) do
                    Literate.markdown(basename(nb); config)
                end
            end
        catch e
            println("An error occured: ", e)
        end
    end
end

any(isnan, ts) && error("Please check errors.")
