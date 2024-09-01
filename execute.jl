using Literate
using JSON
using Pkg
using IJulia
Pkg.activate(Base.current_project())

ENV["GKSwstype"] = "100"

function main(; rmsvg=true)
    file = get(ENV, "NB", "test.ipynb")
    if endswith(file, ".jl")
        run_literate(file; rmsvg)
    elseif endswith(file, ".ipynb")
        IJulia.installkernel("Julia", "--project=@.")
        run_ipynb(file)
    else
        error("$(file) is not a valid notebook file!")
    end
end

# Strip SVG output from a Jupyter notebook
function strip_svg(ipynb)
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

function run_literate(file; rmsvg=true)
    cachedir = get(ENV, "NBCACHE", ".cache")
    outpath = joinpath(abspath(pwd()), cachedir, dirname(file))
    mkpath(outpath)
    ipynb = Literate.notebook(file, outpath; mdstrings=true, execute=true)
    rmsvg && strip_svg(ipynb)
    return
end

function run_ipynb(file)
    cachedir = get(ENV, "NBCACHE", ".cache")
    outpath = joinpath(abspath(pwd()), cachedir, file)
    mkpath(dirname(outpath))
    kernelname = "--ExecutePreprocessor.kernel_name=julia-1.$(VERSION.minor)"
    execute = get(ENV, "ALLOWERRORS", " ") == "true" ? "--execute --allow-errors" : "--execute"
    timeout = "--ExecutePreprocessor.timeout=" * get(ENV, "TIMEOUT", "-1")
    cmd = `jupyter nbconvert --to notebook $(execute) $(timeout) $(kernelname) --output $(outpath) $(file)`
    run(cmd)
    return
end

main()
