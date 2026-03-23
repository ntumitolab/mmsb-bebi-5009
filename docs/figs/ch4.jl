# # Chapter 4
using OrdinaryDiffEq
using Catalyst
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DisplayAs: PNG
using Plots
Plots.gr(linewidth=1.5)

function get_gradient(prob, xsym, ysym, xrange, yrange; t = nothing, normalize=0.9)
    ## The order of state variables (unknowns) in the ODE system is not guaranteed. So we may need to swap the order of x and y when calling ∂F.
    swap_or_not(x, y; xidx=1) = xidx == 1 ? [x, y] : [y, x]
    ∂F = prob.f
    ps = prob.p
    sys = prob.f.sys
    xidx = ModelingToolkit.variable_index(sys, xsym)
    yidx = ModelingToolkit.variable_index(sys, ysym)
    xx = [x for y in yrange, x in xrange]
    yy = [y for y in yrange, x in xrange]
    dx = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[xidx], xx, yy)
    dy = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[yidx], xx, yy)
    if normalize > 0
        maxnormx = maximum(dx)
        maxnormy = maximum(dy)
        dx .*= normalize / maxnormx * step(xrange)
        dy .*= normalize / maxnormy * step(yrange)
    end
    return (; xx, yy, dx, dy)
end

# Steady states and phase plots in an asymmetric network.
@time "Build system" rn402 = @reaction_network begin
    k1 / (1 + B^n), 0 --> A
    k2, 0 --> B
    k3, A --> 0
    k4, B --> 0
    k5, A --> B
end

up402 = Dict(:k1 => 20.0, :k2 => 5.0, :k3 => 5.0, :k4 => 5.0, :k5 => 2.0, :n => 4.0, :A => 0.0, :B => 0.0)


# ## Fig 4.2 A
tend = 1.5
prob402 = ODEProblem(rn402, up402, tend)
u0s = [
    [0.0, 0.0],
    [0.5, 0.6],
    [0.17, 1.1],
    [0.25, 1.9],
    [1.85, 1.7]
]

#---
@time sols = map(u0s) do u0
    solve(remake(prob402, u0=u0), Tsit5())
end
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2 A")

# ### Fig. 4.2 B (Phase plot)
@unpack A, B = rn402
plot(sols[1], idxs=(A, B), xlabel="Time", ylabel="Concentration", title="Fig. 4.2 B", aspect_ratio=1, size=(600, 600), label=false, xlims=(0, 2), ylims=(0, 2))


# ## Fig. 4.3 A (Multiple time series)
plot(xlabel="Time", ylabel="Concentration")
for sol in sols
    plot!(sol, label=false)
end
plot!(title="Fig. 4.3 A")

# ## Fig. 4.3 B (Phase plot)
plot(xlabel="[A]", ylabel="[B]")
for sol in sols
    plot!(sol, idxs=(A, B), label=false)
end
plot!(title="Fig. 4.3 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2))

# Let's sketch vector fields in phase plots
xrange = 0:0.1:2
yrange = 0:0.1:2
(; xx, yy, dx, dy) = get_gradient(prob402, A, B, xrange, yrange)

pl = quiver(xx, yy, quiver=(dx, dy), xlabel="[A]", ylabel="[B]", title="Fig. 4.4 A", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), color=:gray)
for sol in sols
    plot!(pl, sol, idxs=(A, B), label=false)
end
pl

# ## Figure 4.5A
## nullclines
contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline")

## Phase plots
for sol in sols
    plot!(sol, idxs=(A, B), label=false)
end

plot!(title="Fig. 4.5 A", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), xlabel="[A]", ylabel="[B]", legend=:bottomright)


# ## Figure 4.5 B
# Vector field with nullclines
quiver(xx, yy, quiver=(dx, dy), xlabel="[A]", ylabel="[B]", title="Fig. 4.5 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:bottomright)

# ## Fig 4.7 A
# Symmetric (bistable) biological networks.
@time "Build system" rn407 = @reaction_network begin
    k1 / (1 + B^n1), 0 --> A
    k2 / (1 + A^n2), 0 --> B
    k3, A --> 0
    k4, B --> 0
end

up407 = Dict(:k1 => 20.0, :k2 => 20.0, :k3 => 5.0, :k4 => 5.0, :n1 => 4.0, :n2 => 1.0, :A => 3.0, :B => 1.0)

tend = 4.0
@time "Build problem" prob407 = ODEProblem(rn407, up407, (0.0, tend))
@time sol1 = solve(prob407, Tsit5())
@time sol2 = solve(remake(prob407, u0=[:A => 1.0, :B => 3.0]), Tsit5())
pl1 = plot(sol1, xlabel="Time", ylabel="Concentration", title="Fig 4.7A (1)")
pl2 = plot(sol2, xlabel="Time", ylabel="Concentration", title="Fig 4.7A (2)")
plot(pl1, pl2, layout=(2, 1))

# ## Fig 4.7 B
# Vector field with nullclines
xrange = 0:0.01:5
yrange = 0:0.01:5
@unpack A, B = rn407
(; xx, yy, dx, dy) = get_gradient(prob407, A, B, xrange, yrange)

contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false) |> PNG
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false) |> PNG
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:bottomright)
quiver!(xx[1:10:end, 1:10:end], yy[1:10:end, 1:10:end], quiver=(dx[1:10:end, 1:10:end], dy[1:10:end, 1:10:end]), color=:gray)
plot!(title="Fig 4.7 B", aspect_ratio=1, size=(600, 600), xlims=(0, 5), ylims=(0, 5), xlabel="[A]", ylabel="[B]", legend=:bottomright)

# ## Fig 4.8
# Symmetric parameter set
prob408 = remake(prob407, p=[:k1 => 20.0, :k2 => 20.0, :k3 => 5.0, :k4 => 5.0, :n1 => 4.0, :n2 => 4.0], u0=[:A => 3.0, :B => 1.0], tspan=(0.0, 4.0))
@time sol1 = solve(prob408, Tsit5())
@time sol2 = solve(remake(prob408, u0=[:A => 1.0, :B => 3.0]), Tsit5())
pl1 = plot(sol1, xlabel="Time", ylabel="Concentration", title="Fig 4.8A (1)")
pl2 = plot(sol2, xlabel="Time", ylabel="Concentration", title="Fig 4.8A (2)")
plot(pl1, pl2, layout=(2, 1))

# Nullclines and vector field
@unpack A, B = prob408.f.sys
xrange = 0:0.01:5
yrange = 0:0.01:5
(; xx, yy, dx, dy) = get_gradient(prob408, A, B, xrange, yrange)
contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false) |> PNG
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false) |> PNG
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:bottomright)
quiver!(xx[1:10:end, 1:10:end], yy[1:10:end, 1:10:end], quiver=(dx[1:10:end, 1:10:end], dy[1:10:end, 1:10:end]), color=:gray)
plot!(title="Fig 4.8 B", aspect_ratio=1, size=(600, 600), xlims=(0, 5), ylims=(0, 5), xlabel="[A]", ylabel="[B]", legend=:bottomright)

# ## Fig 4.8 C
# Around the unstable steady-state
xrange = range(1.0, 1.5, length=21)
yrange = range(1.0, 1.5, length=21)
(; xx, yy, dx, dy) = get_gradient(prob408, A, B, xrange, yrange)
contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false) |> PNG
plot!([], [],  line=(:black, :solid), label="A nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false) |> PNG
plot!([], [], line=(:black, :dash), label="B nullcline", legend=:bottomright)
quiver!(xx, yy, quiver=(dx, dy), color=:gray)
plot!(title="Fig 4.8 C", aspect_ratio=1, size=(600, 600), xlims=(1.0, 1.5), ylims=(1.0, 1.5), xlabel="[A]", ylabel="[B]", legend=:bottomright)

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

# ## Fig 4.11
# Surface plots
using Plots
z1(x, y) = x^2 + 0.5y^2
z2(x, y) = (.2x^2-1)^2 + y^2
x1 = y1 = range(-1.0, 1.0, 41)
x2 = range(-2.75, 2.75, 41)
y2 = range(-0.75, 0.75, 41)
p1 = surface(x1, y1, z1, title="Single-well potential")
p2 = contourf(x1, y1, z1)
p3 = surface(x2, y2, z2, title="Double-well potential")
p4 = contourf(x2, y2, z2)
plot(p1, p2, p3, p4, size=(800, 600))
