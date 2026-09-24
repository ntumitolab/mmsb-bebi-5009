# # Fig. 7.30
# model of synthetic band detector system
using DifferentialEquations
using SteadyStateDiffEq
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

#---
@time "Build system" rn730 = @reaction_network begin
    hillr(L, aG, KL, 2), 0 --> G
    bG, G --> 0
    hillr(C, aL1, KC, 2), 0 --> L
    mm(R, aL2, KR), 0 --> L
    bL, L --> 0
    mm(R, aC, KR), 0 --> C
    bC, C --> 0
    k1 * (RT - 2R)^2 * A^2, 0 --> R
    k2, R --> 0
end

#---
up730 = Dict(
    :aG => 2.0, # muM/min
    :bG => 0.07, # /min
    :KL => 0.8, # muM
    :aL1 => 1, # muM/min
    :aL2 => 1, # muM/min
    :KC => 0.008, # muM
    :bL => 0.02, # /min
    :aC => 1, # muM/min
    :bC => 0.07, # /min
    :k1 => 0.5, # /muM^3 /min
    :k2 => 0.02, # /min
    :KR => 0.01, # muM
    :RT => 0.5, # muM
    :A => 0.1,
    :G => 0.0,
    :L => 0.0,
    :C => 0.0,
    :R => 0.0,
)

(sols, as) = let N = 1:100
    as = [exp10(-4 + 4 * (i - 1) / length(N)) for i in N]
    @time "Build problem" prob = SteadyStateProblem(rn730, up730)
    prob_func = (prob, ctx) -> remake(prob, p=[:A => as[ctx.sim_id]])
    eprob = EnsembleProblem(prob; prob_func)
    @time "Solve problem" sols = solve(eprob, DynamicSS(FBDF()); trajectories=length(N))
    (sols, as)
end;

luxI = map(s -> s[:L], sols.u)
CI = map(s -> s[:C], sols.u)
GFP = map(s -> s[:G], sols.u)
plot(as, [luxI, CI, GFP], xscale=:log10, xlabel="AHL concentration (μM)", ylabel="Concentration (μM)", title="Fig. 7.30", label=["LuxI" "cI" "GFP"])
