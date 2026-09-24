# # Fig 7.28
# model of synthetic pulse generating system
using DifferentialEquations
using Catalyst
using Plots

function vg728(R, C, KR, KC, aG)
    fR = R / KR
    fC = C / KC
    return aG * (fR / (1 + fR + fC^2 + fR * fC^2))
end

@time "Build system" rn728 = @reaction_network begin
    vg728(R, C, KR, KC, aG), 0 --> G
    bG, G --> 0
    mm(R, aC, KR), 0 --> C
    bC, C --> 0
    k1 * (RT - 2R)^2 * A^2, 0 --> R
    k2, R --> 0
end

#---
up728 = Dict(
    :aG => 80.0, # muM/min
    :bG => 0.07, # /min
    :KC => 0.008, # muM
    :aC => 0.5, # muM/min
    :bC => 0.3, # /min
    :k1 => 0.5, # /muM^3 /min
    :k2 => 0.02, # /min
    :KR => 0.02, # muM
    :RT => 0.5, # muM
    :A => 10.0, # muM
    :G => 0.0,
    :C => 0.0,
    :R => 0.0,
)

alg = FBDF()
tend = 50.0
@time "Build problem" prob728 = ODEProblem(rn728, up728, (0.0, tend))
@time "Solve problem" sol728 = solve(prob728, alg)

plot(sol728, idxs=[:G, :C, :R], xlabel="Time (min)", ylabel="Concentration (μM)", title="Fig 7.28", label=["GFP" "cI" "LuxR:AHL complex"])
