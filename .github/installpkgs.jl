using Pkg, Dates
Pkg.Registry.update()
Pkg.instantiate()
Pkg.precompile()
Pkg.gc(collect_delay=Day(0))
