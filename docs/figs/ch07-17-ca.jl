# # Fig 7.17
# Goodwin oscillator model: https://en.wikipedia.org/wiki/Goodwin_model_(biology)
using ComponentArrays
using OrdinaryDiffEq
using Plots
Plots.gr(linewidth=1.5)

#---
## The reactions, exactly as written in the Catalyst version:
## ```
## (a / (k^n + Z^n), b), 0 <--> X
## (α * X, β),           0 <--> Y
## (γ * Y, δ),           0 <--> Z
## ```
mass_action(k, x) = k * x  ## first-order flux, e.g. first-order removal

## Hill repression of the production of `X` by the end product `Z`: the flux
## falls from `a` to `0` as `Z` grows past the threshold `k`, and `n` sets how
## sharply, i.e. how switch-like the feedback is.
hill_repression(z, p) = p.a / (p.k^p.n + z^p.n)

function goodwin_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1p = hill_repression(u.Z, p)       ## 0 --> X, repressed by Z
    v1m = mass_action(p.b, u.X)         ## X --> 0
    v2p = mass_action(p.α, u.X)         ## 0 --> Y, at a rate proportional to X
    v2m = mass_action(p.β, u.Y)         ## Y --> 0
    v3p = mass_action(p.γ, u.Y)         ## 0 --> Z, at a rate proportional to Y
    v3m = mass_action(p.δ, u.Z)         ## Z --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.X = v1p - v1m
    du.Y = v2p - v2m
    du.Z = v3p - v3m
    return nothing
end

#---
p = ComponentArray(a=360.0, k=1.368, b=1.0, α=1.0, β=0.6, γ=1.0, δ=0.8, n=12.0)
u0 = ComponentArray(X=0.0, Y=0.0, Z=0.0)

alg = FBDF()
tend = 35.0
@time "Build problem" prob717 = ODEProblem(goodwin_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol717 = solve(prob717, alg)

#---
plot(sol717, idxs=[:X, :Y, :Z], labels=["X" "Y" "Z"], title="Fig 7.17 (A)", xlabel="Time", ylabel="Concentration")

#---
## The same trajectory as a phase portrait in the three species.
plot(sol717, idxs=(:X, :Y, :Z), labels=false, title="Fig 7.17 (B)", size=(600, 600), xlabel="X", ylabel="Y", zlabel="Z")
