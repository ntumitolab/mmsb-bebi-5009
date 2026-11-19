# # Fig 6.5
# Model of G-protein signalling pathway
using ComponentArrays
using DiffEqCallbacks
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
using Startup: mass_action, with_param
Plots.gr(linewidth=1.5)

#---
## The reactions
## ```
## kRL * L, R --> RL
## kRLm, RL --> R
## kGa, RL + G --> RL + Ga + Gbg
## kGd0, Ga --> Gd
## kG1, Gd + Gbg --> G
## ```
#---
function gprotein_ode!(du, u, p, t)
    v1 = mass_action(p.kRL, p.L, u.R)     ## R + L --> RL
    v2 = mass_action(p.kRLm, u.RL)        ## RL --> R
    v3 = mass_action(p.kGa, u.RL, u.G)    ## RL + G --> RL + Ga + Gbg
    v4 = mass_action(p.kGd0, u.Ga)        ## Ga --> Gd
    v5 = mass_action(p.kG1, u.Gd, u.Gbg)  ## Gd + Gbg --> G

    du.R = -v1 + v2
    du.RL = v1 - v2
    du.G = -v3 + v5
    du.Ga = v3 - v4
    du.Gd = v4 - v5
    du.Gbg = v3 - v5
    return nothing
end

## The rate constants and the initial condition of the Catalyst model. `L`
## starts OFF, as the discrete parameter `L ~ 0` there.
p = ComponentArray(kRL=2e6, kRLm=0.01, kGa=1e-5, kGd0=0.11, kG1=1.0, L=0.0)
u0 = ComponentArray(R=4e3, RL=0.0, G=1e4, Ga=0.0, Gd=0.0, Gbg=0.0)

# ## Fig 6.5 (A)
L_signal(t) = 1e-9 * (200.0 <= t < 800.0)
update_ligand!(integrator) = (integrator.p.L = L_signal(integrator.t))
cbs = PresetTimeCallback([200.0, 800.0], update_ligand!)

#---
alg = FBDF()
tend = 1200.0
@time "Build problem" prob605 = ODEProblem(gprotein_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol605 = solve(prob605, alg; callback=cbs)

plot(sol605, idxs=[:RL, :Ga], labels=["RL" "Ga"], title="Fig. 6.05 (A)", xlabel="Time", ylabel="Concentration")

# ## Fig 6.5 (B)
lrange = range(0, 20 * 1e-9, length=101)

@time "Build problem" prob605b = SteadyStateProblem(gprotein_ode!, u0, p)
## Sweeping `L` means rebuilding the whole parameter array: pair-based
## `remake(prob; p=[:L => lval])` silently no-ops for a `SteadyStateProblem`.
@time "Ensemble simulations" sols605 = map(lrange) do lval
    _p = remake(prob605b; p=with_param(p, :L, lval))
    sssol = solve(_p, DynamicSS(FBDF()))
end

ga = [sssol.u[:Ga] for sssol in sols605]
rl = [sssol.u[:RL] for sssol in sols605]
plot(lrange .* 1e9, [ga, rl], label=["Ga" "RL"], xlabel="L (nM)", ylabel="Steady-state abundance", title="Fig. 6.5 (B)")
