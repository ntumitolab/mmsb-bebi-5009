# # Fig 6.3 A
# Two component pathway for Figure 6.3 (A).
using ComponentArrays
using DiffEqCallbacks
using OrdinaryDiffEq
using Plots
using Startup: mass_action
Plots.gr(linewidth=1.5)

#---
## The reactions
## ```
## (k1 * L, km1), R <--> RL
## k2, RL + P --> RL + Ps
## k3, Ps --> P
## ```
#---
function pathway_ode!(du, u, p, t)
    v1 = mass_action(p.k1, p.L, u.R)   ## R + L --> RL
    v2 = mass_action(p.km1, u.RL)      ## RL --> R + L
    v3 = mass_action(p.k2, u.RL, u.P)  ## RL + P --> RL + Ps
    v4 = mass_action(p.k3, u.Ps)       ## Ps --> P

    du.R = -v1 + v2
    du.RL = v1 - v2
    du.P = -v3 + v4
    du.Ps = v3 - v4
    return nothing
end

L_signal(t) = 3.0 * (1.0 <= t < 3.0)
update_ligand!(integrator) = (integrator.p.L = L_signal(integrator.t))
cbs = PresetTimeCallback([1.0, 3.0], update_ligand!)

#---
p = ComponentArray(k1=5.0, km1=1.0, k2=6.0, k3=2.0, L=0.0)
u0 = ComponentArray(R=3.0, RL=0.0, P=8.0, Ps=0.0)

alg = FBDF()
tend = 10.0
@time "Build problem" prob603 = ODEProblem(pathway_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol603 = solve(prob603, alg; callback=cbs)

#---
plot(sol603, idxs=[:RL, :Ps], labels=["RL" "P*"], title="Fig. 6.3 (A)", xlabel="Time", ylabel="Concentration")
plot!(sol603.t, L_signal.(sol603.t), st=:steppost, label="L")
