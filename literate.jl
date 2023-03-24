using Literate

cd("docs")

config = Dict("mdstrings" => true)

Literate.notebook("hw-01.jl"; config)
Literate.notebook("intro-01-first-steps.jl"; config)
Literate.notebook("intro-02-plotting.jl"; config)
Literate.notebook("intro-03-diffeq.jl"; config)
Literate.notebook("intro-04-gillespie.jl"; config)
Literate.notebook("mmsb-1.jl"; config)
Literate.notebook("mmsb-2.jl"; config)
Literate.notebook("mmsb-3.jl"; config)
Literate.notebook("mmsb-4.jl"; config)
Literate.notebook("mmsb-5.jl"; config)
