# # Fig 7.23 and 7.31
# Hasty synthetic oscillator model, without Catalyst. The reactions
# ```
# v0_723(x, y, alpha, sigma),            0 --> x
# ay * v0_723(x, y, alpha, sigma),       0 --> y
# gammax,                                x --> 0
# gammay,                                y --> 0
# ```
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: get_gradient, mass_action, meshgrid
Plots.gr(linewidth=1.5)

#---
## The production rate of the oscillator: autocatalytic in `x` through the
## `x^2` and `x^4` terms, and repressed by `y` through `(1 + y^4)`.
v0_723(x, y, alpha, sigma) = (1 + x^2 + alpha * sigma * x^4) / ((1 + x^2 + sigma * x^4) * (1 + y^4))

## Each reaction contributes a *flux* to the balances: the rate that is
## written in front of the reaction, times the mass-action monomial of its
## reactants (`mass_action` is one of the reaction-network helpers in `Startup`).

## Single oscillator: `0 --> x`, `0 --> y`, `x --> 0`, `y --> 0`.
function hasty_ode!(du, u, p, t)
    ## Flux through each reaction:
    v0 = v0_723(u.x, u.y, p.alpha, p.sigma)     ## the shared production rate
    v1 = mass_action(v0)                        ## 0 --> x
    v2 = mass_action(p.ay * v0)                 ## 0 --> y
    v3 = mass_action(p.gammax, u.x)             ## x --> 0
    v4 = mass_action(p.gammay, u.y)             ## y --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.x = v1 - v3
    du.y = v2 - v4
    return nothing
end

# ## Fig. 7.23 (A)
alg = FBDF()
p723 = ComponentArray(alpha=11.0, sigma=2.0, gammax=0.2, gammay=0.012, ay=0.2)
u0723 = ComponentArray(x=0.3963, y=2.3346)

tend = 300.0
@time "Build problem" prob723 = ODEProblem(hasty_ode!, u0723, (0.0, tend), p723)
@time "Solve problem" sol723 = solve(prob723, alg)
plot(sol723, idxs=[:x, :y], xlabel="Time", ylabel="Concentration", title="Fig 7.23 (A)", legend=:right)

# ## Fig. 7.23 (B)
# Vector field
xrange = range(0, 1.5, 51)
yrange = range(0, 3, 51)

xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(hasty_ode!, u0723, p723, :x, :y, xx, yy)
contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [], line=(:black, :solid), label="x nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="y nullcline")
plot!(sol723, idxs=(:x, :y), label="Trajectory")
plot!(xlabel="x", ylabel="y", title="Fig 7.23 (B)", legend=:bottomright, size=(600, 600))

# ## Fig. 7.31
## Two cells, each running the oscillator above, that exchange `x` and `y` by
## diffusion. The pair of reactions `(D, D), x1 <--> x2` is two first-order
## reactions, `D, x1 --> x2` with flux `D * x1` and `D, x2 --> x1` with flux
## `D * x2`, so the net flux from cell 1 to cell 2 is `D * (x1 - x2)`; likewise
## for `y`.
function coupled_hasty_ode!(du, u, p, t)
    ## Flux through each reaction:
    v01 = v0_723(u.x1, u.y1, p.alpha, p.sigma)  ## the production rate in cell 1
    v02 = v0_723(u.x2, u.y2, p.alpha, p.sigma)  ## the production rate in cell 2
    v1 = mass_action(v01)                       ## 0 --> x1
    v2 = mass_action(p.ay * v01)                ## 0 --> y1
    v3 = mass_action(v02)                       ## 0 --> x2
    v4 = mass_action(p.ay * v02)                ## 0 --> y2
    v5 = mass_action(p.gammax, u.x1)            ## x1 --> 0
    v6 = mass_action(p.gammay, u.y1)            ## y1 --> 0
    v7 = mass_action(p.gammax, u.x2)            ## x2 --> 0
    v8 = mass_action(p.gammay, u.y2)            ## y2 --> 0

    ## Net diffusion fluxes out of cell 1, from `x1 <--> x2` and `y1 <--> y2`:
    jx = mass_action(p.D, u.x1) - mass_action(p.D, u.x2)
    jy = mass_action(p.D, u.y1) - mass_action(p.D, u.y2)

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.x1 = v1 - v5 - jx
    du.y1 = v2 - v6 - jy
    du.x2 = v3 - v7 + jx
    du.y2 = v4 - v8 + jy
    return nothing
end

#---
alg = FBDF()
p731 = ComponentArray(alpha=11.0, sigma=2.0, gammax=0.2, gammay=0.012, ay=0.2, D=0.015)
u0731 = ComponentArray(x1=0.3963, y1=2.3346, x2=0.5578, y2=1.9317)

tend = 500.0
@time "Build problem" prob731 = ODEProblem(coupled_hasty_ode!, u0731, (0.0, tend), p731)
@time "Solve problem" sol731 = solve(prob731, alg)

## `plot(sol)` would label the components by their names, but here the four
## curves are listed explicitly so the legend reads `x1, y1, x2, y2`.
plot(sol731, idxs=[:x1, :y1, :x2, :y2], labels=["x1" "y1" "x2" "y2"], xlabel="Time (a.u.)", ylabel="Concentration (a.u.)", title="Fig. 7.31", legend=:left)
