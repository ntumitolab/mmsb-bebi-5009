# # Fig 7.25
# Model of quorum sensing mechanism of Vibrio fischeri.
using ComponentArrays
using DifferentialEquations
using SteadyStateDiffEq
using Plots
using Startup: hil, with_param
Plots.gr(linewidth=1.5)

#---
## Hill functions: `hil(x, k)` is the saturating `x / (k + x)` and the
## three-argument form adds a Hill coefficient `n`. Both come from `Startup`.

# Model with/without feedback
function make_ode(feedback::Bool)
    function qs_ode!(du, u, p, t)
        R0 = p.RT - 2u.Rstar                     ## free LuxR, algebraic
        v0 = feedback ? p.k0 * u.I : p.k0 * 15.0 ## 0 --> A, made by LuxI
        v1 = p.k1 * u.A^2 * R0^2                 ## 2A + 2R0 --> Rstar
        v2 = p.k2 * u.Rstar                      ## Rstar --> 2A + 2R0
        va = p.n * (u.A - u.Aout)                ## A --> Aout by diffusion

        du.A = v0 - 2v1 + 2v2 - va
        du.I = p.a0 + p.a * hil(u.Rstar, p.KM) - p.b * u.I
        du.Rstar = v1 - v2
        du.Aout = p.popsize * va - p.diff * u.Aout
        return nothing
    end
    return qs_ode!
end

#---
p = ComponentArray(
    k0 = 0.0008, k1 = 0.5, k2 = 0.02, n = 0.6, a = 10.0, b = 0.07,
    a0 = 0.05, KM = 0.01, RT = 0.5, diff = 1000.0, popsize = 1000.0,
)
u0 = ComponentArray(A = 0.0, I = 0.0, Rstar = 0.0, Aout = 0.0)

@time "Build problem" prob725 = SteadyStateProblem(make_ode(true), u0, p)
@time "Build problem" prob725_no_feedback = SteadyStateProblem(make_ode(false), u0, p)

#---
## Sweeping `popsize` means rebuilding the whole parameter array: pair-based
## `remake(prob; p=[:popsize => n])` silently no-ops for a `SteadyStateProblem`.
npops = 1:50:5001
trajectories = length(npops)
alg = DynamicSS(FBDF())
prob_func = (prob, ctx) -> remake(prob; p = with_param(prob.p, :popsize, npops[ctx.sim_id]))
eprob = EnsembleProblem(prob725; prob_func)
eprob_nofeed = EnsembleProblem(prob725_no_feedback; prob_func)
@time sim = solve(eprob, alg; trajectories);
@time sim_no_feed = solve(eprob_nofeed, alg; trajectories);

luxI = map(s -> s.u[:I], sim.u)
luxI_nofeed = map(s -> s.u[:I], sim_no_feed.u)
plot(npops, [luxI, luxI_nofeed], xscale=:log10, xlabel="Population size", ylabel="LuxI concentration (μM)", title="Fig. 7.25", label=["With feedback" "Without feedback"], legend=:topleft)
