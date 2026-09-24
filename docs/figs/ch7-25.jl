# # Fig 7.25
# Model of quorum sensing mechanism of Vibrio fischeri
using DifferentialEquations
using SteadyStateDiffEq
using Catalyst
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

# Model with/without feedback
function build_sys725(; feedback=true, simplify=true, name=:model725)
    hil(x, k) = x / (k + x)
    hil(x, k, n) = hil(x^n, k^n)
    @parameters k0 = 0.0008 k1 = 0.5 k2 = 0.02 n = 0.6 a = 10 b = 0.07 a0 = 0.05 KM = 0.01 RT = 0.5 diff = 1000 popsize = 1000
    @variables A(t) = 0 I(t) = 0 Rstar(t) = 0 Aout(t) = 0 R0(t)
    R0 = RT - 2Rstar
    v0 = feedback ? k0 * I : k0 * 15.0
    v1 = k1 * A^2 * R0^2
    v2 = k2 * Rstar
    va = n * (A - Aout)
    eqs = [
        R0 ~ RT - 2Rstar,
        D(A) ~ v0 - 2v1 + 2v2 - va,
        D(I) ~ a0 + a * hil(Rstar, KM) - b * I,
        D(Rstar) ~ v1 - v2,
        D(Aout) ~ popsize * va - diff * Aout
    ]
    sys = ODESystem(eqs, t; name)
    return simplify ? mtkcompile(sys) : sys
end

#---

@time "Build system" sys725 = build_sys725()
@time "Build system" sys725_no_feedback = build_sys725(feedback=false)
@time "Build problem" prob725 = SteadyStateProblem(sys725, [])
@time "Build problem" prob725_no_feedback = SteadyStateProblem(sys725_no_feedback, [])

npops = 1:50:5001
trajectories = length(npops)
alg = DynamicSS(FBDF())
@unpack popsize, I = sys725
prob_func = (prob, ctx) -> remake(prob, p=[popsize => npops[ctx.sim_id]])
eprob = EnsembleProblem(prob725; prob_func)
eprob_nofeed = EnsembleProblem(prob725_no_feedback; prob_func)
@time sim = solve(eprob, alg; trajectories);
@time sim_no_feed = solve(eprob_nofeed, alg; trajectories);

luxI = map(s -> s[I], sim.u)
luxI_nofeed = map(s -> s[I], sim_no_feed.u)
plot(npops, [luxI, luxI_nofeed], xscale=:log10, xlabel="Population size", ylabel="LuxI concentration (μM)", title="Fig. 7.25", label=["With feedback" "Without feedback"], legend=:topleft)
