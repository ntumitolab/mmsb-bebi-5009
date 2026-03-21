# # Chapter 7
# ## Fig 7.25
# model of quorum sensing mechanism of Vibrio fischeri
using OrdinaryDiffEq
using SteadyStateDiffEq
using Catalyst
using Plots

# Model with feedback
@time "Build system" rn725 = @reaction_network begin
    @combinatoric_ratelaws false
    k0 * I, 0 --> A
    (k1, k2), 2A + 2R0 <--> Rstar
    n * (A - Aout), A => 0
    (a0 + mm(Rstar, a, KM), b), 0 <--> I
    popsize * n * (A - Aout), 0 --> Aout
    diffAout, Aout --> 0
end

# Model without feedback
@time "Build system" rn725_no_feedback = @reaction_network begin
    @combinatoric_ratelaws false
    15k0, 0 --> A
    (k1, k2), 2A + 2R0 <--> Rstar
    n * (A - Aout), A => 0
    (a0 + mm(Rstar, a, KM), b), 0 <--> I
    popsize * n * (A - Aout), 0 --> Aout
    diffAout, Aout --> 0
end

up725 = Dict(
    :k0=>0.0008, # /min
    :k1=>0.5, # /muM^3 /min
    :k2=>0.02, # /min
    :n=>0.6,
    :a=>10, # muM/min
    :b=>0.07, # /min
    :a0=>0.05, # muM/min
    :KM=>0.01, # muM
    :R0=>0.5, # muM
    :diffAout=>1000,
    :popsize=>1000,
    :A=>0.0,
    :I=>0.0,
    :Rstar=>0.0,
    :Aout=>0.0,
)

@time "Build problem" prob725 = SteadyStateProblem(rn725, up725)
@time "Build problem" prob725_no_feedback = SteadyStateProblem(rn725_no_feedback, up725)

# npops = 1:50:5000
# trajectories = length(npops)
# alg = DynamicSS(KenCarp47())
# prob_func = (prob, i, repeat) -> remake(prob, p=[:popsize=>npops[i]])
# eprob = EnsembleProblem(prob725; prob_func)
# eprob_nofeed = EnsembleProblem(prob725_no_feedback; prob_func)
# @time sim = solve(eprob, alg; trajectories);
# @time sim_nofeed = solve(eprob_nofeed, alg; trajectories);



# ## Fig 7.28
# model of synthetic pulse generating system
using OrdinaryDiffEq
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

up728 = Dict(
    :aG=>80.0, # muM/min
    :bG=>0.07, # /min
    :KC=>0.008, # muM
    :aC=>0.5, # muM/min
    :bC=>0.3, # /min
    :k1=>0.5, # /muM^3 /min
    :k2=>0.02, # /min
    :KR=>0.02, # muM
    :RT=>0.5, # muM
    :A=>10.0, # muM
    :G=>0.0,
    :C=>0.0,
    :R=>0.0,
)

tend = 50.0
@time "Build problem" prob728 = ODEProblem(rn728, up728, (0.0, tend))
@time "Solve problem" sol728 = solve(prob728, KenCarp47())

# Fig 7.28
plot(sol728, idxs=[:G, :C, :R], xlabel="Time (min)", ylabel="Concentration (μM)", title="Fig 7.28", label=["GFP" "cI" "LuxR:AHL complex"])

# ## Fig. 7.30
# model of synthetic band detector system
using OrdinaryDiffEq
using SteadyStateDiffEq
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

# Model
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

up730 = Dict(
    :aG=>2.0, # muM/min
    :bG=>0.07, # /min
    :KL=>0.8, # muM
    :aL1=>1, # muM/min
    :aL2=>1, # muM/min
    :KC=>0.008, # muM
    :bL=>0.02, # /min
    :aC=>1, # muM/min
    :bC=>0.07, # /min
    :k1=>0.5, # /muM^3 /min
    :k2=>0.02, # /min
    :KR=>0.01, # muM
    :RT=>0.5, # muM
    :A=>0.1,
    :G=>0.0,
    :L=>0.0,
    :C=>0.0,
    :R=>0.0,
)

(sols, as) = let N= 1:100
    as = [exp10(-4 + 4*(i-1)/length(N)) for i in N]
    @time "Build problem" prob = SteadyStateProblem(rn730, up730)
    prob_func = (prob, i, repeat) -> remake(prob, p=[:A=>as[i]])
    eprob = EnsembleProblem(prob; prob_func)
    @time "Solve problem" sols = solve(eprob, DynamicSS(KenCarp47()); trajectories=length(N))
    (sols, as)
end

# Figure
luxI = map(s -> s[:L], sols)
CI = map(s -> s[:C], sols)
GFP = map(s -> s[:G], sols)
plot(as, [luxI, CI, GFP], xscale=:log10, xlabel="AHL concentration (μM)", ylabel="Concentration (μM)", title="Fig. 7.30", label=["LuxI" "cI" "GFP"])


# ## Fig. 7.31
# Hasty synthetic oscillator model
using OrdinaryDiffEq
using Catalyst
using Plots
Plots.gr(linewidth=1.5)
# Model
v0_731(x, y, alpha, sigma) = (1 + x^2 + alpha * sigma * x^4) / ((1 + x^2 + sigma * x^4) * (1 + y^4))

@time "Build system" rn731 = @reaction_network begin
    v0_731(x1, y1, alpha, sigma), 0 --> x1
    ay * v0_731(x1, y1, alpha, sigma), 0 --> y1
    v0_731(x2, y2, alpha, sigma), 0 --> x2
    ay * v0_731(x2, y2, alpha, sigma), 0 --> y2
    gammax, x1 --> 0
    gammay, y1 --> 0
    gammax, x2 --> 0
    gammay, y2 --> 0
    (D, D), x1 <--> x2
    (D, D), y1 <--> y2
end

up731 = Dict(
    :alpha=>11.0,
    :sigma=>2.0,
    :gammax=>0.2,
    :gammay=>0.012,
    :ay=>0.2,
    :D=>0.015,
    :x1=>0.3963,
    :y1=>2.3346,
    :x2=>0.5578,
    :y2=>1.9317,
)

tend = 500.0
@time "Build problem" prob731 = ODEProblem(rn731, up731, (0.0, tend))
@time "Solve problem" sol731 = solve(prob731, KenCarp47())

# Fig 7.31
plot(sol731, xlabel="Time (a.u.)", ylabel="Concentration (a.u.)", title="Fig. 7.31", legend=:left)

# ## Fig. 7.38
# stochastic implementation (Gillespie SSA) of a collection of spontaneously decaying molecules.
using JumpProcesses
using Catalyst
using Plots

@time "Build system" rn738 = @reaction_network begin
    d, A --> 0
end

tend = 5.0
ps738 = Dict(:d=>1.0)
@time "Build problem" jprob1000 = JumpProblem(rn738, [:A=>1000], (0.0, tend), ps738)
jprob100 = remake(jprob1000, u0=[:A=>100])
jprob10 = remake(jprob1000, u0=[:A=>10])

@time "Solve problem" sol1000 = solve(jprob1000)
@time "Solve problem" sol100 = solve(jprob100)
@time "Solve problem" sol10 = solve(jprob10)

# Fig 7.38
plot(sol1000, xlabel="Time", ylabel="# of molecules", title="Fig. 7.38 (A)", label="1000 molecules")
plot!(t->1000 * exp(-t), linestyle=:dash, label="ODE solution")

#---
plot(sol100, xlabel="Time", ylabel="# of molecules", title="Fig. 7.38 (B)", label="100 molecules")
plot!(t->100 * exp(-t), linestyle=:dash, label="ODE solution")

#---
plot(sol10, xlabel="Time", ylabel="# of molecules", title="Fig. 7.38 (C)", label="10 molecules")
plot!(t->10 * exp(-t), linestyle=:dash, label="ODE solution")
