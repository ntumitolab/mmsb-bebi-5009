using Distributed
using PrettyTables
using SHA
using IJulia

@everywhere begin
    ENV["GKSwstype"] = "100"
    using Literate, Pkg, JSON
end

# Strip SVG output from a Jupyter notebook
@everywhere function strip_svg(ipynb)
    @info "Stripping SVG in $(ipynb)"
    nb = open(JSON.parse, ipynb, "r")
    for cell in nb["cells"]
        !haskey(cell, "outputs") && continue
        for output in cell["outputs"]
            !haskey(output, "data") && continue
            datadict = output["data"]
            if haskey(datadict, "image/png") || haskey(datadict, "image/jpeg")
                delete!(datadict, "text/html")
                delete!(datadict, "image/svg+xml")
            end
        end
    end
    rm(ipynb)
    open(ipynb, "w") do io
        JSON.print(io, nb, 1)
    end
    return ipynb
end

# Remove cached notebook and sha files if there is no corresponding notebook
function clean_cache(cachedir)
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
end

function list_notebooks(basedir, cachedir)
    ipynbs = String[]
    litnbs = String[]

    for (root, dirs, files) in walkdir(basedir)
        for file in files
            if endswith(file, ".ipynb") || endswith(file, ".jl")
                nb = joinpath(root, file)
                shaval = read(nb, String) |> sha256 |> bytes2hex
                @info "Notebook $(nb): hash=$(shaval)"
                shafilename = joinpath(cachedir, root, splitext(file)[1] * ".sha")
                if isfile(shafilename) && read(shafilename, String) == shaval
                    @info "Notebook $(nb) cache hits and will not be executed."
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
    return (; ipynbs, litnbs)
end

@everywhere function run_literate(file, cachedir; rmsvg=true)
    outpath = joinpath(abspath(pwd()), cachedir, dirname(file))
    mkpath(outpath)
    ipynb = Literate.notebook(file, outpath; mdstrings=true, execute=true)
    rmsvg && strip_svg(ipynb)
    return ipynb
end

function main(;
    basedir=get(ENV, "DOCDIR", "docs"),
    cachedir=get(ENV, "NBCACHE", ".cache"),
    printtable=true, rmsvg=true)

    mkpath(cachedir)
    clean_cache(cachedir)

    (; ipynbs, litnbs) = list_notebooks(basedir, cachedir)

    # Execute literate notebooks in worker process(es)
    ts_lit = pmap(litnbs; on_error=ex -> NaN) do nb
        @elapsed run_literate(nb, cachedir; rmsvg)
    end
    rmprocs(workers()) # Remove worker processes to release some memory

    # Debug notebooks one by one if there are errors
    for (nb, t) in zip(litnbs, ts_lit)
        if isnan(t)
            println("Debugging notebook: ", nb)
            try
                withenv("JULIA_DEBUG" => "Literate") do
                    run_literate(nb, cachedir; rmsvg)
                end
            catch e
                println(e)
            end
        end
    end

    any(isnan, ts_lit) && error("Please check literate notebook error(s).")

    # Install IJulia kernel
    IJulia.installkernel("Julia", "--project=@.", "--heap-size-hint=3G")

    # nbconvert command array
    ntasks = parse(Int, get(ENV, "NBCONVERT_JOBS", "1"))
    kernelname = "--ExecutePreprocessor.kernel_name=julia-1.$(VERSION.minor)"
    execute = ifelse(get(ENV, "ALLOWERRORS", " ") == "true", "--execute --allow-errors", "--execute")
    timeout = "--ExecutePreprocessor.timeout=" * get(ENV, "TIMEOUT", "-1")
    cmds = [`jupyter nbconvert --to notebook $(execute) $(timeout) $(kernelname) --output $(joinpath(abspath(pwd()), cachedir, nb)) $(nb)` for nb in ipynbs]

    # Run the nbconvert commands in parallel
    ts_ipynb = asyncmap(cmds; ntasks) do cmd
        @elapsed run(cmd)
    end

    # Print execution result
    printtable && pretty_table([litnbs ts_lit; ipynbs ts_ipynb], header=["Notebook", "Elapsed (s)"])
end

# Run code
main()
