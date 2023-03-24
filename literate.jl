using Literate

config = Dict("mdstrings" => true)

Literate.notebook("docs/hw-01.jl", "docs/"; config)
Literate.notebook("docs/intro-01-first-steps.jl", "docs/"; config)
Literate.notebook("docs/intro-02-plotting.jl", "docs/"; config)
Literate.notebook("docs/intro-03-diffeq.jl", "docs/"; config)
Literate.notebook("docs/intro-04-gillespie.jl", "docs/"; config)
Literate.notebook("docs/mmsb-1.jl", "docs/"; config)
Literate.notebook("docs/mmsb-2.jl", "docs/"; config)
Literate.notebook("docs/mmsb-3.jl", "docs/"; config)
Literate.notebook("docs/mmsb-4.jl", "docs/"; config)
Literate.notebook("docs/mmsb-5.jl", "docs/"; config)
