using Pkg, Dates
Pkg.add(["PrettyTables", "Literate"])
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
Pkg.gc(collect_delay=Day(0))
