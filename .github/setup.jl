using Pkg
using Dates
Pkg.Registry.update()
Pkg.add(["JSON", "Literate"])
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
if ENV["RUNNER_ENVIRONMENT"] == "github-hosted"
    Pkg.gc(;collect_delay=Day(0))
end
using Literate
