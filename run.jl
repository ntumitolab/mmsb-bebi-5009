
using Distributed
using PrettyTables
using SHA
using IJulia

@everywhere begin
    ENV["GKSwstype"] = "100"
    using Literate, Pkg
    Pkg.activate(Base.current_project())
end

outfile = "ipynbs.txt"
basedir = get(ENV, "DOCDIR", "docs") # Defaults to docs/
cachedir = get(ENV, "NBCACHE", ".cache") # Defaults to .cache/
ipynbs = String[]
litnbs = String[]

# Collect the list of notebooks
for (root, dirs, files) in walkdir(basedir)
    for file in files
        if endswith(file, ".ipynb") || endswith(file, ".jl")
            nb = joinpath(root, file)
            shaval = read(nb, String) |> sha256 |> bytes2hex
            @info "Notebook $(nb): hash=$(shaval)"
            shafilename = joinpath(cachedir, root, splitext(file)[1] * ".sha")
            # Cache hit
            if isfile(shafilename) && read(shafilename, String) == shaval
                @info "Notebook $(nb) cache hits and will not be executed."
            # Cache miss
            else
                @info "Notebook $(nb) cache misses. Writing hash to $(shafilename)."
                mkpath(dirname(shafilename))
                write(shafilename, shaval)
                if endswith(file, ".ipynb")
                    push!(ipynbs, nb)
                elseif endswith(file, ".jl")
                    push!(litnbs, nb)
                end
            end
        end
    end
end

# Remove cached notebook and sha files if there is no respective notebook
for (root, dirs, files) in walkdir(cachedir)
    for file in files
        if endswith(file, ".ipynb") || endswith(file, ".sha")
            fn = joinpath(joinpath(splitpath(root)[2:end]), splitext(file)[1])
            nb = fn * ".ipynb"
            lit = fn * ".jl"
            if !isfile(nb) && !isfile(lit)
                fullfn = joinpath(root, file)
                @info "Notebook $(nb) or $(lit) not found. Removing $(fullfn)."
                rm(fullfn)
            end
        end
    end
end

# Execute literate notebooks in worker process(es)
ts = pmap(litnbs; on_error=ex->NaN) do nb
    outdir = joinpath(cachedir, dirname(nb))
    @elapsed Literate.notebook(nb, outdir; mdstrings=true)
end

# Show literate notebook execution results
pretty_table([litnbs ts], header=["Notebook", "Elapsed (s)"])

# Remove worker processes in Distributed.jl
rmprocs(workers())

# Debug notebooks one by one if there are errors
for (nb, t) in zip(litnbs, ts)
    if isnan(t)
        println("Debugging notebook: ", nb)
        try
            withenv("JULIA_DEBUG" => "Literate") do
                Literate.notebook(nb, dirname(nb); mdstrings=true)
            end
        catch e
            println(e)
        end
    end
end

any(isnan, ts) && error("Please check literate notebook error(s).")

# Install IJulia kernel
IJulia.installkernel("Julia", "--project=@.", "--heap-size-hint=3G")

# nbconvert command array
ntasks = parse(Int, get(ENV, "NBCONVERT_JOBS", "1"))
kernelname = "--ExecutePreprocessor.kernel_name=julia-1.$(VERSION.minor)"
execute = ifelse(get(ENV, "ALLOWERRORS", " ") == "true", "--execute --allow-errors", "--execute")
timeout = "--ExecutePreprocessor.timeout=" * get(ENV, "TIMEOUT", "-1")
cmds = [`jupyter nbconvert --to notebook $(execute) $(timeout) $(kernelname) --output $(joinpath(abspath(pwd()), cachedir, nb)) $(nb)` for nb in ipynbs]

# Run the nbconvert commands in parallel
ts = asyncmap(cmds; ntasks) do cmd
    @elapsed run(cmd)
end

# Print execution result
pretty_table([ipynbs ts], header=["Notebook", "Elapsed (s)"])
