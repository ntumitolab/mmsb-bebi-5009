# # Fig 7.21
# Repressilator model: three genes whose proteins repress each other in a
# cycle, `A` represses the mRNA of `B`, `B` that of `C`, and `C` that of `A`.
using ComponentArrays
using OrdinaryDiffEq
using Plots
Plots.gr(linewidth=1.5)

#---
mass_action(k, x) = k * x  ## first-order flux, e.g. first-order removal

## Hill repression `hillr(x, α, 1, n)` of Catalyst: the production rate falls
## from `α` to `0` as the repressor `x` grows past the threshold `K = 1`, and
## `n` sets how sharply, i.e. how switch-like the repression is. The general
## form is
## ```
## α * K^n / (x^n + K^n)
## ```
## which with `K = 1` reduces to `α / (1 + x^n)`.
hill_repression(x, p) = p.α * 1.0^p.n / (x^p.n + 1.0^p.n)

## Production of an mRNA: a small basal rate `α0` plus the repressed maximal
## rate, so the gene is never completely off.
mrna_production(x, p) = p.α0 + hill_repression(x, p)

function repressilator_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1p = mrna_production(u.C, p)      ## 0 --> mA, repressed by C
    v1m = mass_action(1.0, u.mA)       ## mA --> 0
    v2p = mrna_production(u.A, p)      ## 0 --> mB, repressed by A
    v2m = mass_action(1.0, u.mB)       ## mB --> 0
    v3p = mrna_production(u.B, p)      ## 0 --> mC, repressed by B
    v3m = mass_action(1.0, u.mC)       ## mC --> 0
    v4p = mass_action(p.β, u.mA)       ## 0 --> A, at a rate proportional to mA
    v4m = mass_action(p.β, u.A)        ## A --> 0
    v5p = mass_action(p.β, u.mB)       ## 0 --> B, at a rate proportional to mB
    v5m = mass_action(p.β, u.B)        ## B --> 0
    v6p = mass_action(p.β, u.mC)       ## 0 --> C, at a rate proportional to mC
    v6m = mass_action(p.β, u.C)        ## C --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.mA = v1p - v1m
    du.mB = v2p - v2m
    du.mC = v3p - v3m
    du.A = v4p - v4m
    du.B = v5p - v5m
    du.C = v6p - v6m
    return nothing
end

#---
p = ComponentArray(α=298.2, α0=0.03, n=2.0, β=0.2)
u0 = ComponentArray(mA=0.2, mB=0.3, mC=0.4, A=0.1, B=0.1, C=0.5)

alg = FBDF()
tend = 300.0
@time "Build problem" prob721 = ODEProblem(repressilator_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol721 = solve(prob721, alg)

#---
plot(sol721, idxs=[:A, :B, :C], xlabel="Time", ylabel="Concentration", title="Fig 7.21", legend=:left)
