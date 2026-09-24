# # Fig 6.19
# Sine wave response of g-protein signalling pathway
using OrdinaryDiffEq
using Catalyst
using Plots
using LaTeXStrings

@time "Build system" rn605 = @reaction_network begin
    @parameters lt l_AMP l_per
    @equations begin
        L ~ lt + (lt / l_AMP) * cospi(2t / l_per)
    end
    kRL * L, R --> RL
    kRLm, RL --> R
    kGa, RL + G --> RL + Ga + Gbg
    kGd0, Ga --> Gd
    kG1, Gd + Gbg --> G
end

#---
alg = FBDF()
up619 = Dict(
    :kRL => 2e6,
    :kRLm => 0.01,
    :kGa => 1e-5,
    :kGd0 => 0.11,
    :kG1 => 1.0,
    :R => 4e3,
    :G => 1e4,
    :lt => 1e-9,
    :l_per => 200.0,
    :l_AMP => 5.0,
    :RL => 0.0,
    :Ga => 0.0,
    :Gd => 0.0,
    :Gbg => 0.0
)

tend = 1000.0
@time "Build problem" prob619 = ODEProblem(rn605, up619, (0.0, tend); mtkcompile=true)
@time "Solve problem" sol619 = solve(prob619, alg)

@unpack Ga, L = rn605
plot(sol619, idxs=[Ga, L*1E12], title="Fig 6.19 (A)", xlabel="Time", ylabel="Abundance", label=["Ga" L"L \cdot 10^{12}"])
