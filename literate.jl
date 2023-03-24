using Literate

cd("docs")

config = Dict("mdstrings" => true)

notebooks = (
    "hw-01.jl",
    "intro-01-first-steps.jl",
    "intro-02-plotting.jl",
    "intro-03-diffeq.jl",
    "intro-04-gillespie.jl",
    "mmsb-1.jl",
    "mmsb-2.jl",
    "mmsb-3.jl",
    "mmsb-4.jl",
    "mmsb-5.jl"
)

for nb in notebooks
    Literate.notebook(nb; config)
end
