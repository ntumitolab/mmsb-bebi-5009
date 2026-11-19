# # Fig 2.09
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: mass_action
Plots.gr(linewidth=1.5)

#---
## The reactions
## ```
## k1, 0 --> A
## k2, A --> B
## k3, A + B --> C + D
## k4, C --> 0
## k5, D --> 0
## ```
function metabolic_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1 = mass_action(p.k1)              ## 0 --> A
    v2 = mass_action(p.k2, u.A)         ## A --> B
    v3 = mass_action(p.k3, u.A, u.B)    ## A + B --> C + D
    v4 = mass_action(p.k4, u.C)         ## C --> 0
    v5 = mass_action(p.k5, u.D)         ## D --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.A = v1 - v2 - v3
    du.B = v2 - v3
    du.C = v3 - v4
    du.D = v3 - v5
    return nothing
end

#---
p = ComponentArray(k1=3.0, k2=2.0, k3=2.5, k4=3.0, k5=4.0)
u0 = ComponentArray(A=0.0, B=0.0, C=0.0, D=0.0)

tend = 10.0
@time "Build problem" prob = ODEProblem(metabolic_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol = solve(prob, Tsit5())

#---
## Component names work as indices into the solution too: `idxs=[:A, :B, :C, :D]`
## selects the species by name, without knowing their order in `u`.
plot(sol, idxs=[:A, :B, :C, :D], labels=["A" "B" "C" "D"], title="Fig 2.9", xlabel="Time", ylabel="Concentration")
