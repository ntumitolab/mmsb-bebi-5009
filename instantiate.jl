using Pkg, Dates

Pkg.add(["IJulia", "Literate", "PrettyTables"])
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
Pkg.gc(collect_delay=Day(0))
