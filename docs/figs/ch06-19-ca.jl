# # Fig 6.19 (A)
# Sine wave response of g-protein signalling pathway
using ComponentArrays
using OrdinaryDiffEq
using Plots
using LaTeXStrings
Plots.gr(linewidth=1.5)

#---
mass_action(k) = k                ## zero-order flux
mass_action(k, x) = k * x         ## first-order flux
mass_action(k, x, y) = k * x * y  ## second-order flux

L_signal(t, p) = p.lt + (p.lt / p.l_AMP) * cospi(2t / p.l_per)

#---
function gprotein_ode!(du, u, p, t)
    L = L_signal(t, p)            ## the ligand, sampled from its definition

    v1 = mass_action(p.kRL, L, u.R)       ## R + L --> RL
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

#---
## The rate constants and the initial condition of the Catalyst model.
p = ComponentArray(kRL=2e6, kRLm=0.01, kGa=1e-5, kGd0=0.11, kG1=1.0, lt=1e-9, l_AMP=5.0, l_per=200.0)
u0 = ComponentArray(R=4e3, RL=0.0, G=1e4, Ga=0.0, Gd=0.0, Gbg=0.0)

alg = FBDF()
tend = 1000.0
@time "Build problem" prob619 = ODEProblem(gprotein_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol619 = solve(prob619, alg)

#---
plot(sol619, idxs=[:Ga], title="Fig 6.19 (A)", xlabel="Time", ylabel="Abundance", label="Ga")
plot!(sol619.t, 1e12 .* L_signal.(sol619.t, Ref(p)), label=L"L \cdot 10^{12}")
