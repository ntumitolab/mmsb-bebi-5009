# # Fig 4.7
# Symmetric (bistable) biological networks.
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: meshgrid, normalize_gradient
Plots.gr(linewidth=1.5)

#---
## The reactions are
## ```
## k1 / (1 + B^n1), 0 --> A
## k2 / (1 + A^n2), 0 --> B
## k3, A --> 0
## k4, B --> 0
## ```
function toggle_ode!(du, u, p, t)
    du.A = p.k1 / (1 + u.B^p.n1) - p.k3 * u.A
    du.B = p.k2 / (1 + u.A^p.n2) - p.k4 * u.B
    return nothing
end

function get_gradient(f!, p, xx, yy)
    dx = map(xx, yy) do x, y
        u = ComponentArray(A = x, B = y)
        du = similar(u)
        f!(du, u, p, 0.0)
        du.A
    end
    dy = map(xx, yy) do x, y
        u = ComponentArray(A = x, B = y)
        du = similar(u)
        f!(du, u, p, 0.0)
        du.B
    end
    return dx, dy
end

# ## Fig 4.7 A
alg = FBDF()
p407 = ComponentArray(k1=20.0, k2=20.0, k3=5.0, k4=5.0, n1=4.0, n2=1.0)
u407 = ComponentArray(A=3.0, B=1.0)
tend = 4.0
@time "Build problem" prob407 = ODEProblem(toggle_ode!, u407, (0.0, tend), p407)
@time sol1 = solve(prob407, alg)
@time sol2 = solve(remake(prob407, u0=[:A => 1.0, :B => 3.0]), alg)
pl1 = plot(sol1, idxs=[:A, :B], labels=["A" "B"], xlabel="Time", ylabel="Concentration", title="Fig 4.7A (1)")
pl2 = plot(sol2, idxs=[:A, :B], labels=["A" "B"], xlabel="Time", ylabel="Concentration", title="Fig 4.7A (2)")
plot(pl1, pl2, layout=(2, 1))

# ## Fig 4.7 B
# Vector field with nullclines
xrange = range(0, 5, 201)
yrange = range(0, 5, 201)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(toggle_ode!, prob407.p, xx, yy)
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
plot!(title="Fig. 4.7B", aspect_ratio=1, size=(600, 600), xlims=(0, 5), ylims=(0, 5), xlabel="[A]", ylabel="[B]", legend=:topright)

# ## Fig 4.8 A
# Symmetric parameter set
prob408 = remake(prob407, p=[:k1 => 20.0, :k2 => 20.0, :k3 => 5.0, :k4 => 5.0, :n1 => 4.0, :n2 => 4.0], u0=[:A => 3.0, :B => 1.0], tspan=(0.0, 4.0))
@time sol1 = solve(prob408, alg)
@time sol2 = solve(remake(prob408, u0=[:A => 1.0, :B => 3.0]), alg)
pl1 = plot(sol1, idxs=[:A, :B], labels=["A" "B"], xlabel="Time", ylabel="Concentration", title="Fig 4.8A (1)")
pl2 = plot(sol2, idxs=[:A, :B], labels=["A" "B"], xlabel="Time", ylabel="Concentration", title="Fig 4.8A (2)")
plot(pl1, pl2, layout=(2, 1))

# ## Fig 4.8 B
# Nullclines and vector field
xrange = range(0, 5, 201)
yrange = range(0, 5, 201)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(toggle_ode!, prob408.p, xx, yy)
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
plot!(title="Fig. 4.8 B", aspect_ratio=1, size=(600, 600), xlims=(0, 5), ylims=(0, 5), xlabel="[A]", ylabel="[B]", legend=:topright)

# ## Fig 4.8 C
# Around the unstable steady-state.
xrange = 1.0:0.01:1.5
yrange = 1.0:0.01:1.5
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(toggle_ode!, prob408.p, xx, yy)
xx_sp = xx[1:5:end, 1:5:end]
yy_sp = yy[1:5:end, 1:5:end]
dx_sp = dx[1:5:end, 1:5:end]
dy_sp = dy[1:5:end, 1:5:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange)*5, step(yrange)*5)
quiver(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline")
plot!(title="Fig. 4.8 C", aspect_ratio=1, size=(600, 600), xlims=(1.0, 1.5), ylims=(1.0, 1.5), xlabel="[A]", ylabel="[B]", legend=:topright)

# Another way to draw nullclines is to find the analytical solutions for dA (or dB) is zero. And then sketch the nullclines in a parameteric plot.
nca47(b, p) = p.k1 / p.k3 / (1 + b^p.n1)
ncb47(a, p) = p.k2 / p.k4 / (1 + a^p.n2)

pls = map((8.0, 16.0, 20.0, 35.0)) do k1
    ps = (k1=k1, k2=20., k3=5., k4=5., n1=4., n2=4.)
    pl = plot(t->nca47(t, ps), identity, 0, 7, color=:red, label="Nullcline A")
    plot!(pl, identity, t->ncb47(t, ps), 0, 7, color=:blue, label="Nullcline B")
    plot!(xlims=(0, 7), ylims=(0, 7), aspect_ratio=1, title="k1 = $k1", xlabel="[A]", ylabel="[B]", legend=:right)
end

plot(pls..., layout=(2, 2), size=(800, 800))
