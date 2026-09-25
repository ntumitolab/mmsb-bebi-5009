# # Fig 4.15
# Oscillatory networks
using Startup: meshgrid, get_gradient, normalize_gradient
using Catalyst
using DiffEqCallbacks
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
Plots.gr(linewidth=1.5)

# ## Fig 4.15 A
@time "Build system" model415 = @reaction_network begin
    k0 , 0 --> A
    k1 * (1 + B^n), A --> B
    k2, B --> 0
end

#---
alg = FBDF()
up415 = Dict(:k0 => 8.0, :k1 => 1.0, :k2 => 5.0, :n => 2.0, :A => 1.5, :B => 1.0)
tend = 8.0
@time "Build problem" prob415 = ODEProblem(model415, up415, (0.0, tend))

u0s = [
    [:A=>1.5, :B=>1.0],
    [:A=>0.0, :B=>1.0],
    [:A=>0.0, :B=>3.0],
    [:A=>2.0, :B=>0.0],
]

@time "Solve problems" sols = map(u0s) do u0
    solve(remake(prob415, u0=u0), alg)
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.15 A")

# ## Fig 4.15 B
# Vector field with nullclines
xrange = range(0, 4, 101)
yrange = range(0, 4, 101)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(prob415, :A, :B, xx, yy)
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
dx, dy = get_gradient(prob416, :A, :B, xx, yy)
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
prob417 = remake(prob416, u0=[:A=>2.0, :B=>1.5], tspan=(0.0, 10.0))
@time sol = solve(prob417, alg)

plot(sol, idxs=(:A, :B), arrow=(:head))

xrange = range(1, 3, 201)
yrange = range(1, 3, 201)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(prob417, :A, :B, xx, yy)
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
