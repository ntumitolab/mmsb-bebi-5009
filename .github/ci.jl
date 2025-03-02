using Distributed
using Tables
using MarkdownTables
using SHA
using IJulia

@everywhere begin
    ENV["GKSwstype"] = "100"
    using Literate, Pkg, JSON
end

# Strip SVG output from a Jupyter notebook
@everywhere function strip_svg(ipynb)
    oldfilesize = filesize(ipynb)
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
    write(ipynb, JSON.json(nb, 1))
    @info "Stripped SVG in $(ipynb). The original size is $(oldfilesize). The new size is $(filesize(ipynb))."
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

"Recursively list Jupyter and Literate notebooks. Also process caching."
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

# Run a Literate.jl notebook
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
    rmsvg=true)

    mkpath(cachedir)
    clean_cache(cachedir)

    (; ipynbs, litnbs) = list_notebooks(basedir, cachedir)

    if !isempty(litnbs)
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
    else
        ts_lit = []
    end

    if !isempty(ipynbs)
        IJulia.installkernel("Julia", "--project=@.")

        # nbconvert command array
        ntasks = parse(Int, get(ENV, "NBCONVERT_JOBS", "1"))
        kernelname = "--ExecutePreprocessor.kernel_name=julia-1.$(VERSION.minor)"
        execute = ifelse(get(ENV, "ALLOWERRORS", " ") == "true", "--execute --allow-errors", "--execute")
        timeout = "--ExecutePreprocessor.timeout=" * get(ENV, "TIMEOUT", "-1")

        # Run the nbconvert commands in parallel
        ts_ipynb = asyncmap(ipynbs; ntasks) do nb
            @elapsed begin
                nbout = joinpath(abspath(pwd()), cachedir, nb)
                cmd = `jupyter nbconvert --to notebook $(execute) $(timeout) $(kernelname) --output $(nbout) $(nb)`
                run(cmd)
                rmsvg && strip_svg(nbout)
            end
        end
    else
        ts_ipynb = []
    end

    # Print execution result
    Tables.table([litnbs ts_lit; ipynbs ts_ipynb]; header=["Notebook", "Elapsed (s)"]) |> markdown_table(String) |> print
end

# Run code
main()
