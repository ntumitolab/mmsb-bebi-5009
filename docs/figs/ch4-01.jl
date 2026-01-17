#===
# Fig 4.2, 4.3, 4.4, 4.5

Steady states and phase plots in an asymmetric network.
===#
using OrdinaryDiffEq
import ComponentArrays as CA
using SimpleUnPack
using CairoMakie

# The model for figure 4.1, 4.2, and 4.3.
_dA401(A, B, p, t) = p.k1 / (1 + B^p.n) - (p.k3 + p.k5) * A
_dB401(A, B, p, t) = p.k2 + p.k5 * A - p.k4 * B

function model401!(D, u, p, t)
    @unpack A, B = u
    D.A = _dA401(A, B, p, t)
    D.B = _dB401(A, B, p, t)
    return nothing
end

# ## Fig 4.2 A
tend = 1.5
ps402a = CA.ComponentArray(
    k1 = 20.0,
    k2 = 5.0,
    k3 = 5.0,
    k4 = 5.0,
    k5 = 2.0,
    n = 4.0
)
ics402a = CA.ComponentArray(
    A = 0.0,
    B = 0.0
)
prob401 = ODEProblem(model401!, ics402a, tend, ps402a)

u0s = [
    CA.ComponentArray(
        A = a,
        B = b
    ) for (a, b) in [
    [0.0, 0.0],
    [0.5, 0.6],
    [0.17, 1.1],
    [0.25, 1.9],
    [1.85, 1.7]]
]

#---
@time sols = map(u0s) do u0
    solve(remake(prob401, u0=u0), Tsit5())
end

#---
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig. 4.2 A\nTime series"
)
lines!(ax, 0..tend, t-> sols[1](t).A, label = "A")
lines!(ax, 0..tend, t-> sols[1](t).B, label = "B")
axislegend(ax, position = :rt)
fig

# ## Fig. 4.2 B (Phase plot)
fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.2 B\nPhase plot",
    aspect = 1,
)
let ts = 0:0.01:tend
    aa = [sols[1](t).A for t in ts]
    bb = [sols[1](t).B for t in ts]
    lines!(ax, aa, bb)
end
limits!(ax, 0.0, 2.0, 0.0, 2.0)
fig

# ## Fig. 4.3 A (Multiple time series)
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig. 4.3A\nMultiple time series"
)
for sol in sols
    lines!(ax, 0..tend, t-> sol(t).A, alpha=0.7)
    lines!(ax, 0..tend, t-> sol(t).B, alpha=0.7)
end
fig

# ## Fig. 4.3 B (Phase plot)
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.3 B\nPhase plot",
    aspect = 1,
)
for sol in sols
    ts = 0:0.01:tend
    aa = [sol(t).A for t in ts]
    bb = [sol(t).B for t in ts]
    lines!(ax, aa, bb)
end
limits!(ax, 0.0, 2.0, 0.0, 2.0)
fig

# Let's sketch vector fields in phase plots
∂F44 = function (x, y)
    da = _dA401(x, y, ps402a, nothing)
    db = _dB401(x, y, ps402a, nothing)
    Point2d(da, db)
end

fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.4 A\nPhase plot with vector field",
    aspect = 1,
)

streamplot!(ax, ∂F44, 0..2, 0..2)
for sol in sols
    ts = 0:0.01:tend
    aa = [sol(t).A for t in ts]
    bb = [sol(t).B for t in ts]
    lines!(ax, aa, bb, color=:black)
end
limits!(ax, 0.0, 2.0, 0.0, 2.0)
fig

# ## Figure 4.5A
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.5 A\nPhase plot with nullclines",
    aspect = 1,
)

## Phase plots
for sol in sols
    ts = 0:0.01:tend
    aa = [sol(t).A for t in ts]
    bb = [sol(t).B for t in ts]
    lines!(ax, aa, bb, color=:black)
end

## nullclines
xs = 0:0.01:2
ys = 0:0.01:2
zA44 = [_dA401(x, y, ps402a, nothing) for x in xs, y in ys]
zB44 = [_dB401(x, y, ps402a, nothing) for x in xs, y in ys]
contour!(ax, xs, ys, zA44, levels=[0], color=:red)
lines!(ax, Float64[], Float64[], color=:red, label="A nullcline")
contour!(ax, xs, ys, zB44, levels=[0], color=:blue)
lines!(ax, Float64[], Float64[], color=:blue, label="B nullcline")
limits!(ax, 0.0, 2.0, 0.0, 2.0)
axislegend(ax, position = :rb)
fig

# ## Figure 4.5 B
# Vector field
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.5 B\nVector field with nullclines",
    aspect = 1,
)
## Nullclines
contour!(ax, xs, ys, zA44, levels=[0], color=:red)
lines!(ax, Float64[], Float64[], color=:red, label="A nullcline")
contour!(ax, xs, ys, zB44, levels=[0], color=:blue)
lines!(ax, Float64[], Float64[], color=:blue, label="B nullcline")
streamplot!(ax, ∂F44, 0..2, 0..2)
limits!(ax, 0.0, 2.0, 0.0, 2.0)
fig
