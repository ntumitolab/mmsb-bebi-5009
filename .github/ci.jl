using Distributed
using Tables
using MarkdownTables
using SHA

@everywhere begin
    ENV["GKSwstype"] = "100"
    using Literate, JSON
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
    rm(nbpath; force=true)
    write(nbpath, JSON.json(nb, 1))
    @info "Stripped SVG in $(nbpath). The original size is $(Base.format_bytes(oldfilesize)). The new size is $(Base.format_bytes(filesize(nbpath)))."
    return nbpath
end

# Remove cached notebook and sha files if there is no corresponding notebook
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

# Convert a Jupyter notebook into a Literate notebook. Adapted from https://github.com/JuliaInterop/NBInclude.jl.
function to_literate(nbpath; shell_or_help=r"^\s*[;?]")
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

# List notebooks without caches in a file tree
function list_notebooks(basedir, cachedir)
    list = String[]
    for (root, _, files) in walkdir(basedir)
        for file in files
            name, ext = splitext(file)
            if ext == ".ipynb" || ext == ".jl"
                nb = joinpath(root, file)
                shaval = read(nb, String) |> sha256 |> bytes2hex
                @info "$(nb) SHA256 = $(shaval)"
                shafilename = joinpath(cachedir, root, name * ".sha")
                if isfile(shafilename) && read(shafilename, String) == shaval
                    @info "$(nb) cache hits and will not be executed."
                else
                    @info "$(nb) cache misses. Writing hash to $(shafilename)."
                    mkpath(dirname(shafilename))
                    write(shafilename, shaval)
                    if ext == ".ipynb"
                        litnb = to_literate(nb)
                        rm(nb; force=true)
                        push!(list, litnb)
                    elseif ext == ".jl"
                        push!(list, nb)
                    end
                end
            end
        end
    end
    return list
end

# Run a Literate notebook
@everywhere function run_literate(file, cachedir; rmsvg=true)
    outpath = joinpath(abspath(pwd()), cachedir, dirname(file))
    mkpath(outpath)
    ipynb = Literate.notebook(file, dirname(file); mdstrings=true, execute=true)
    rmsvg && strip_svg(ipynb)
    cp(ipynb, joinpath(outpath, basename(ipynb)); force=true)
    return ipynb
end

function main(;
    basedir=get(ENV, "DOCDIR", "docs"),
    cachedir=get(ENV, "NBCACHE", ".cache"),
    rmsvg=true)

    mkpath(cachedir)
    clean_cache(cachedir)
    litnbs = list_notebooks(basedir, cachedir)

    if !isempty(litnbs)
        # Execute literate notebooks in worker process(es)
        ts_lit = pmap(litnbs; on_error=identity) do nb
            @elapsed run_literate(nb, cachedir; rmsvg)
        end
        rmprocs(workers()) # Remove worker processes to release some memory
        failed = false
        for (nb, t) in zip(litnbs, ts_lit)
            if t isa ErrorException
                println("Notebook: ", nb, "failed with error: \n", t.msg)
                failed = true
            end
        end

        if failed
            error("Please check literate notebook error(s).")
        else
            # Print execution result
            Tables.table([litnbs ts_lit]; header=["Notebook", "Elapsed (s)"]) |> markdown_table(String) |> print
        end
    end
end

# Run code
main()
