# # Fig 4.2 A
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: get_gradient, mass_action, meshgrid, normalize_gradient
Plots.gr(linewidth=1.5)

# Bistable autoregulatory network for Figures 4.2 - 4.5, without Catalyst.
# The reactions
# ```
# k1 / (1 + B^n), 0 --> A
# k2, 0 --> B
# k3, A --> 0
# k4, B --> 0
# k5, A --> B
# ```
# are translated into ODEs by hand: `B` represses its own activator `A`
# through a Hill term, so `A` has two stable steady states and the network
# is a switch. Both the state and the parameter vector are `ComponentArray`s,
# so species and rate constants are labelled by name (`u.A`, `p.k1`) instead
# of by numeric index.

#---
function bistable_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1 = mass_action(p.k1 / (1 + u.B^p.n))  ## 0 --> A, inhibited by B
    v2 = mass_action(p.k2)                  ## 0 --> B
    v3 = mass_action(p.k3, u.A)             ## A --> 0
    v4 = mass_action(p.k4, u.B)             ## B --> 0
    v5 = mass_action(p.k5, u.A)             ## A --> B

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.A = v1 - v3 - v5
    du.B = v2 - v4 + v5
    return nothing
end

#---
p = ComponentArray(k1=20.0, k2=5.0, k3=5.0, k4=5.0, k5=2.0, n=4.0)
u0 = ComponentArray(A=0.0, B=0.0)

tend = 1.5
@time "Build problem" prob402 = ODEProblem(bistable_ode!, u0, (0.0, tend), p)

## Five starting points, chosen around the two stable steady states: some of
## them fall on one side of the separatrix, some on the other.
u0s = [
    ComponentArray(A=0.0, B=0.0),
    ComponentArray(A=0.5, B=0.6),
    ComponentArray(A=0.17, B=1.1),
    ComponentArray(A=0.25, B=1.9),
    ComponentArray(A=1.85, B=1.7)
]

@time sols = map(u0s) do u0
    solve(remake(prob402, u0=u0), FBDF())
end

## `plot(sol)` draws every component, labelled by its name.
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2 A")

# ## Fig 4.2 B
# Phase plot. `idxs=(:A, :B)` picks the components by name.
plot(sols[1], idxs=(:A, :B), xlabel="Time", ylabel="Concentration", title="Fig. 4.2 B", aspect_ratio=1, size=(600, 600), label=false, xlims=(0, 2), ylims=(0, 2))

# ## Fig 4.3 A
# Multiple time series
plot(xlabel="Time", ylabel="Concentration")
for sol in sols
    plot!(sol, label=false)
end
plot!(title="Fig. 4.3 A")

# ## Fig 4.3 B
# Phase plot
plot(xlabel="[A]", ylabel="[B]")
for sol in sols
    plot!(sol, idxs=(:A, :B), label=false)
end
plot!(title="Fig. 4.3 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2))

# ## Fig 4.4 A
# Let's sketch vector fields in phase plots.
xrange = 0:0.1:2
yrange = 0:0.1:2

xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(bistable_ode!, u0, p, :A, :B, xx, yy)
dx_norm, dy_norm = normalize_gradient(dx, dy, step(xrange), step(yrange))
quiver(xx, yy, quiver=(dx_norm, dy_norm), color=:gray)

for sol in sols
    plot!(sol, idxs=(:A, :B), label=false)
end

plot!(title="Fig. 4.4 A", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2))

# ## Fig 4.5A
# Nullclines
xrange = 0:0.01:2
yrange = 0:0.01:2
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(bistable_ode!, u0, p, :A, :B, xx, yy)

contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline")
for sol in sols
    plot!(sol, idxs=(:A, :B), label=false)
end

plot!(title="Fig. 4.5 A", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), xlabel="[A]", ylabel="[B]", legend=:bottomright)

# ## Fig 4.5 B
# Vector field with nullclines.
xrange = 0:0.01:2
yrange = 0:0.01:2
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(bistable_ode!, u0, p, :A, :B, xx, yy)

xx_sp = xx[1:10:end, 1:10:end]
yy_sp = yy[1:10:end, 1:10:end]
dx_sp = dx[1:10:end, 1:10:end]
dy_sp = dy[1:10:end, 1:10:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange)*10, step(yrange)*10)
quiver(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline")
plot!(title="Fig. 4.5 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), xlabel="[A]", ylabel="[B]", legend=:bottomright)
