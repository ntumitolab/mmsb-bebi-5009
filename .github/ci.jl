using Distributed

@everywhere begin
    ENV["GKSwstype"] = "100"
    using Literate, JSON, SHA
end

# Strip SVG output from a Jupyter notebook
@everywhere function strip_svg(nbpath)
    oldfilesize = filesize(nbpath)
    nb = open(JSON.parse, nbpath, "r")
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
    write(nbpath, JSON.json(nb, 1))
    @info "Stripped SVG in $(nbpath). The original size is $(oldfilesize). The new size is $(filesize(nbpath))."
    return nbpath
end

"Remove cached notebook and sha files if there is no corresponding notebook"
function clean_cache(cachedir)
    for (root, _, files) in walkdir(cachedir)
        for file in files
            fn, ext = splitext(file)
            if ext == ".sha"
                target = joinpath(joinpath(splitpath(root)[2:end]), fn)
                nb = target * ".ipynb"
                lit = target * ".jl"
                if !isfile(nb) && !isfile(lit)
                    cachepath = joinpath(root, fn)
                    @info "Notebook $(nb) or $(lit) not found. Removing $(cachepath) SHA and notebook."
                    rm(cachepath * ".sha")
                    rm(cachepath * ".ipynb"; force=true)
                end
            end
        end
    end
end

"Convert a Jupyter notebook into a Literate notebook. Adapted from https://github.com/JuliaInterop/NBInclude.jl."
function to_literate(nbpath; shell_or_help = r"^\s*[;?]")
    nb = open(JSON.parse, nbpath, "r")
    jlpath = splitext(nbpath)[1] * ".jl"
    open(jlpath, "w") do io
        separator = ""
        for cell in nb["cells"]
            if cell["cell_type"] == "code"
                s = join(cell["source"])
                isempty(strip(s)) && continue # Jupyter doesn't number empty cells
                occursin(shell_or_help, s) && continue  # Skip cells with shell and help commands
                print(io, separator, "#---\n", s)  # Literate code block mark
                separator = "\n\n"
            elseif cell["cell_type"] == "markdown"
                txt = join(cell["source"])
                print(io, separator, "#===\n", txt, "\n===#")
                separator = "\n\n"
            end
        end
    end
    return jlpath
end

"Recursively list Literate notebooks"
function list_notebooks(basedir)
    litnbs = String[]
    for (root, _, files) in walkdir(basedir)
        for file in files
            nb = joinpath(root, file)
            name, ext = splitext(file)
            if ext == ".jl"
                push!(litnbs, nb)
            elseif ext == ".ipynb"
                lit = to_literate(nb)
                rm(nb)
                push!(litnbs, lit)
            end
        end
    end
    return litnbs
end

# Run a Literate notebook
@everywhere function run_literate(file, cachedir; rmsvg=true)
    shaval = read(file, String) |> sha256 |> bytes2hex
    @info "$(file) SHA256 = $(shaval)"
    shafilename = joinpath(cachedir, splitext(file)[1] * ".sha")
    ipynb = joinpath(cachedir, splitext(file)[1] * ".ipynb")
    if isfile(shafilename) && read(shafilename, String) == shaval && isfile(ipynb)
        @info "$(file) cache hits. The notebooks is $(ipynb). It will not be executed."
        return ipynb
    end
    @info "$(file) cache misses. Writing hash to $(shafilename)."
    mkpath(dirname(shafilename))
    write(shafilename, shaval)
    outpath = joinpath(abspath(pwd()), cachedir, dirname(file))
    mkpath(outpath)
    @time "$(file) took" ipynb = Literate.notebook(file, outpath; mdstrings=true, execute=true)
    return rmsvg ? strip_svg(ipynb) : ipynb
end

function main(;
    basedir=get(ENV, "DOCDIR", "docs"),
    cachedir=get(ENV, "NBCACHE", ".cache"),
    rmsvg=true)

    mkpath(cachedir)
    clean_cache(cachedir)
    litnbs = list_notebooks(basedir)

    # Execute literate notebooks in worker process(es)
    ts_lit = pmap(litnbs; on_error=ex -> NaN) do nb
        @elapsed run_literate(nb, cachedir; rmsvg)
    end
    # Remove worker processes to release some memory
    rmprocs(workers())
    # Debug notebooks one by one if there are errors
    for (nb, t) in zip(litnbs, ts_lit)
        if isnan(t)
            println("Debugging notebook: ", nb)
            try
                withenv("JULIA_DEBUG" => "Literate") do
                    Literate.notebook(nb; mdstrings=true, execute=true)
                end
            catch e
                println(e)
            end
        end
    end
    any(isnan, ts_lit) && error("Please check notebook error(s).")
end

# Run code
main()
