# # Fig 7.07
# Model of lac operon in E. coli.
using ComponentArrays
using DiffEqCallbacks
using DifferentialEquations
using SteadyStateDiffEq
using Plots
using Startup: hil, with_param
Plots.gr(linewidth=1.5)

#---
function lac_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1p = p.a1 / (1 + p.RToverK1 * (p.K2 / (p.K2 + u.L))^4)   ## 0 --> M
    v1m = p.δM * u.M                                          ## M --> 0
    v2p = p.c1 * u.M                                          ## 0 --> Y
    v2m = p.δY * u.Y                                          ## Y --> 0
    v3 = (p.kL * u.Y) * hil(p.Le, p.KML)                      ## 0 --> L, saturating in the external lactose Le
    v4 = p.δL * u.L                                           ## L --> 0
    v5 = (u.Y / 4) * (2 * p.kg) * hil(u.L, p.KMg)             ## L => 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.M = v1p - v1m
    du.Y = v2p - v2m
    du.L = v3 - v4 - v5
    return nothing
end

cbs = let
    event_le1 = PresetTimeCallback([500.0], (integrator) -> integrator.p.Le = 50)
    event_le2 = PresetTimeCallback([1000.0], (integrator) -> integrator.p.Le = 100)
    event_le3 = PresetTimeCallback([1500.0], (integrator) -> integrator.p.Le = 150)
    event_le4 = PresetTimeCallback([2000.0], (integrator) -> integrator.p.Le = 0)
    CallbackSet(event_le1, event_le2, event_le3, event_le4)
end

#---
p = ComponentArray(
    δM = 0.48, δY = 0.03, δL = 0.02,
    a1 = 0.29, K2 = 2.92e6, RToverK1 = 213.2, c1 = 18.8,
    kL = 6e4, KML = 680.0, kg = 3.6e3, KMg = 7e5,
    Le = 0.0,
)
u0 = ComponentArray(M = 0.01, Y = 0.1, L = 0.0)

alg = FBDF()

# ## Fig 7.07 (A)
tend = 2500.0
@time "Build problem" prob707a = ODEProblem(lac_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol = solve(prob707a, alg, callback=cbs)

plot(sol, idxs=[:Y], labels="β-galactosidase monomer", title="Fig 7.07 (A)", xlabel="Time (min)", ylabel="Concentration")
plot!(t -> 50 * (500 <= t < 1000) + 100 * (1000 <= t < 1500) + 150 * (1500 <= t < 2000), 0, tend, linestyle=:dash, label="External lactose")

# ## Fig 7.07 (B)
# Comparing the original model and the modified model
## The modified model (`rn707b`) keeps the `M` and `Y` equations but lets a
## fixed amount of enzyme `Enz` produce and consume the inducer `L`:
## ```
## Enz * mm(Le, 4kL, KML), 0 --> L
## δL,                    L --> 0
## Enz * mm(L, 2kg, KMg),  L => 0
## ```
## so neither flux depends on `Y` any more; in particular the consumption
## keeps no `Y / 4` factor. The extra parameter `Enz` is added to a copy of
## `p` below, giving `p_mod`.
function lac_ode_modified!(du, u, p, t)
    ## Flux through each reaction:
    v1p = p.a1 / (1 + p.RToverK1 * (p.K2 / (p.K2 + u.L))^4)   ## 0 --> M
    v1m = p.δM * u.M                                          ## M --> 0
    v2p = p.c1 * u.M                                          ## 0 --> Y
    v2m = p.δY * u.Y                                          ## Y --> 0
    v3 = p.Enz * (4 * p.kL) * hil(p.Le, p.KML)                ## 0 --> L
    v4 = p.δL * u.L                                           ## L --> 0
    v5 = p.Enz * (2 * p.kg) * hil(u.L, p.KMg)                 ## L => 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.M = v1p - v1m
    du.Y = v2p - v2m
    du.L = v3 - v4 - v5
    return nothing
end

#---
## This variant adds a new parameter, so the whole array is rebuilt;
## `with_param` only replaces a parameter that `p` already has.
p_mod = ComponentArray(merge(NamedTuple(p), (Enz = 40.0,)))

## Sweeping `Le` means rebuilding the whole parameter array: pair-based
## `remake(prob; p=[:Le => le])` silently no-ops for a `SteadyStateProblem`.
prob707a = SteadyStateProblem(lac_ode!, u0, p)
prob707b = SteadyStateProblem(lac_ode_modified!, u0, p_mod)
lerange = 0:100
prob_func = (prob, ctx) -> remake(prob; p = with_param(prob.p, :Le, lerange[ctx.sim_id]))
output_func = (sol, ctx) -> (sol.u[:Y] / 4, false)
ssalg = DynamicSS(FBDF())

eprob = EnsembleProblem(prob707a; prob_func, output_func)
eprob_mod = EnsembleProblem(prob707b; prob_func, output_func)
@time sim = solve(eprob, ssalg; trajectories = length(lerange))
@time sim_mod = solve(eprob_mod, ssalg; trajectories = length(lerange))

plot(lerange, [sim.u sim_mod.u], labels=["Original model" "Modified model"], title="Fig 7.07 (B)", xlabel="External lactose concentration", ylabel="β-galactosidase tetramer concentration")
