using Pkg
using Dates
Pkg.Registry.update()
Pkg.add(["JSON", "Literate", "MarkdownTables", "SHA", "Tables"])
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
if ENV["RUNNER_ENVIRONMENT"] == "github-hosted"
    Pkg.gc(;collect_delay=Day(0))
end
using Literate
try
    using PythonPlot
catch e
    @info "PythonPlot not installed. Skipping matplotlib installation."
end
