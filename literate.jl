using Literate

cd("docs")

@time Literate.notebook("hw-01.jl")
@time Literate.notebook("intro-01-first-steps.jl")
@time Literate.notebook("intro-02-plotting.jl")
@time Literate.notebook("intro-03-diffeq.jl")
@time Literate.notebook("intro-04-gillespie.jl")
@time Literate.notebook("mmsb-1.jl")
@time Literate.notebook("mmsb-2.jl")
@time Literate.notebook("mmsb-3.jl")
@time Literate.notebook("mmsb-4.jl")
@time Literate.notebook("mmsb-5.jl")
