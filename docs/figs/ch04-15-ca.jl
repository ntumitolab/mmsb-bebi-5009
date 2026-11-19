# # Fig 4.15
# Oscillatory networks for Figures 4.15 - 4.17, without Catalyst.
# The reactions
# ```
# k0,           0 --> A
# k1 * (1 + B^n), A --> B
# k2,             B --> 0
# ```
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: get_gradient, mass_action, meshgrid, normalize_gradient
Plots.gr(linewidth=1.5)

#---
function oscillatory_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1 = mass_action(p.k0)                        ## 0 --> A
    v2 = mass_action(p.k1 * (1 + u.B^p.n), u.A)   ## A --> B, activated by B
    v3 = mass_action(p.k2, u.B)                   ## B --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.A = v1 - v2
    du.B = v2 - v3
    return nothing
end

#---
# ## Fig 4.15 A
alg = FBDF()
p = ComponentArray(k0=8.0, k1=1.0, k2=5.0, n=2.0)
u0 = ComponentArray(A=1.5, B=1.0)
tend = 8.0
@time "Build problem" prob415 = ODEProblem(oscillatory_ode!, u0, (0.0, tend), p)

## Four starting points, spread around the oscillation:
u0s = [
    ComponentArray(A=1.5, B=1.0),
    ComponentArray(A=0.0, B=1.0),
    ComponentArray(A=0.0, B=3.0),
    ComponentArray(A=2.0, B=0.0),
]

@time "Solve problems" sols = map(u0s) do u0
    solve(remake(prob415, u0=u0), alg)
end

## `plot(sol)` draws every component, labelled by its name.
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.15 A")

# ## Fig 4.15 B
# Vector field with nullclines
xrange = range(0, 4, 101)
yrange = range(0, 4, 101)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(oscillatory_ode!, u0, p, :A, :B, xx, yy)
xx_sp = xx[1:5:end, 1:5:end]
yy_sp = yy[1:5:end, 1:5:end]
dx_sp = dx[1:5:end, 1:5:end]
dy_sp = dy[1:5:end, 1:5:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange)*5, step(yrange)*5)
quiver(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)

for sol in sols
    plot!(sol, idxs=(:A, :B), label=false)
end

contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:bottomright)
plot!(title="Fig. 4.15 B", aspect_ratio=1, size=(600, 600), xlims=(0, 4), ylims=(0, 4), xlabel="[A]", ylabel="[B]", legend=:bottomright)

# ## Fig 4.16 A
# Oscillatory parameter set
prob416 = remake(prob415, p=[:n=>2.5], tspan=(0.0, 100.0))
@time "Solve problems" sols = map(u0s) do u0
    solve(remake(prob416, u0=u0), alg)
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.16 A", tspan=(0.0, 10.0))

# ## Fig 4.16 B
xrange = range(0, 4, 201)
yrange = range(0, 4, 201)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(oscillatory_ode!, prob416.u0, prob416.p, :A, :B, xx, yy)
xx_sp = xx[1:10:end, 1:10:end]
yy_sp = yy[1:10:end, 1:10:end]
dx_sp = dx[1:10:end, 1:10:end]
dy_sp = dy[1:10:end, 1:10:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange)*10, step(yrange)*10)
quiver(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)

for sol in sols
    plot!(sol, idxs=(:A, :B), label=false)
end

contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:bottomright)
plot!(title="Fig. 4.16 B", aspect_ratio=1, size=(600, 600), xlims=(0, 4), ylims=(0, 4), xlabel="[A]", ylabel="[B]", legend=:topright)

# ## Fig 4.17
prob417 = remake(prob416, u0=ComponentArray(A=2.0, B=1.5), tspan=(0.0, 10.0))
@time sol = solve(prob417, alg)

plot(sol, idxs=(:A, :B), arrow=(:head))

xrange = range(1, 3, 201)
yrange = range(1, 3, 201)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(oscillatory_ode!, prob417.u0, prob417.p, :A, :B, xx, yy)
xx_sp = xx[1:10:end, 1:10:end]
yy_sp = yy[1:10:end, 1:10:end]
dx_sp = dx[1:10:end, 1:10:end]
dy_sp = dy[1:10:end, 1:10:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange)*10, step(yrange)*10)
quiver!(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:topright)
plot!(title="Fig. 4.16 C", aspect_ratio=1, size=(600, 600), xlims=(1, 3), ylims=(1, 3), xlabel="[A]", ylabel="[B]", legend=:right)
