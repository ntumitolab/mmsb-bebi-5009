# # Fig 4.15, 4.16, 4.17
# Oscillatory networks.
# ## Figure 4.15 (A)
using OrdinaryDiffEq
import ComponentArrays as CA
using SimpleUnPack
using CairoMakie

# The model
_dA415(A, B, p, t) = p.k0 - p.k1 * A * (1 + B^p.n)
_dB415(A, B, p, t) = p.k1 * A * (1 + B^p.n) - p.k2 * B
function model415!(D, u, p, t)
    @unpack A, B = u
    D.A = _dA415(A, B, p, t)
    D.B = _dB415(A, B, p, t)
    nothing
end

#---
ps415 = CA.ComponentArray(
    k0 = 8.0,
    k1 = 1.0,
    k2 = 5.0,
    n = 2.0
)
u0415 = CA.ComponentArray(
    A = 1.5,
    B = 1.0
)

tend = 8.0
prob415 = ODEProblem(model415!, u0415, (0.0, tend), ps415)

#---
u0s = [
    CA.ComponentArray(A=1.5, B=1.0),
    CA.ComponentArray(A=0.0, B=1.0),
    CA.ComponentArray(A=0.0, B=3.0),
    CA.ComponentArray(A=2.0, B=0.0),
]

@time sols = map(u0s) do u0
    solve(remake(prob415, u0=u0), Tsit5())
end

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig. 4.15 A\nTime series"
)
lines!(ax, 0..tend, t-> sols[1](t).A, label = "A")
lines!(ax, 0..tend, t-> sols[1](t).B, label = "B")
axislegend(ax, position = :rt)
fig

# ## Fig 4.15 (B)
∂F415  = function (x, y)
    da = _dA415(x, y, ps415, nothing)
    db = _dB415(x, y, ps415, nothing)
    return Point2d(da, db)
end

## Grid points
xx = 0:0.01:4
yy = 0:0.01:4
ts = 0:0.05:tend
∂A415 = [_dA415(x, y, ps415, nothing) for x in xx, y in yy]
∂B415 = [_dB415(x, y, ps415, nothing) for x in xx, y in yy]

fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.15 B\nPhase plot",
    aspect = 1,
)

streamplot!(ax, ∂F415, xx, yy)
contour!(ax, xx, yy, ∂A415, levels=[0], color=:red)
lines!(ax, Float64[], Float64[], color=:red, label="A nullcline")
contour!(ax, xx, yy, ∂B415, levels=[0], color=:blue)
lines!(ax, Float64[], Float64[], color=:blue, label="B nullcline")
for sol in sols
    lines!(ax, sol, idxs=(1, 2), color=:tomato, label=nothing)
end
axislegend(ax, position = :rc)
limits!(ax, 0.0, 4.0, 0.0, 4.0)
fig

#===
## Fig 4.16 A

Oscillatory parameter set
===#
ps416 = CA.ComponentArray(
    k0 = 8.0,
    k1 = 1.0,
    k2 = 5.0,
    n = 2.5
)

tend = 100.0
u0s = [
    CA.ComponentArray(A=1.5, B=1.0),
    CA.ComponentArray(A=0.0, B=1.0),
    CA.ComponentArray(A=0.0, B=3.0),
    CA.ComponentArray(A=2.0, B=0.0),
]

prob416 = ODEProblem(model415!, u0415, tend, ps416)

#---
@time sols = map(u0s) do u0
    solve(remake(prob416, u0=u0), Tsit5())
end

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig. 4.16 A\nTime series"
)
lines!(ax, 0..8, t-> sols[1](t).A, label = "A")
lines!(ax, 0..8, t-> sols[1](t).B, label = "B")
axislegend(ax, position = :rt)
fig

# ## Fig 4.16 b
∂F416 = function (x, y)
    da = _dA415(x, y, ps416, nothing)
    db = _dB415(x, y, ps416, nothing)
    return Point2d(da, db)
end

xx = 0:0.01:4
yy = 0:0.01:4

∂A416 = [_dA415(x, y, ps416, nothing) for x in xx, y in yy]
∂B416 = [_dB415(x, y, ps416, nothing) for x in xx, y in yy];
#---
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.16 B\nPhase plot",
    aspect = 1,
)

streamplot!(ax, ∂F416, xx, yy)
contour!(ax, xx, yy, ∂A416, levels=[0], color=:red)
lines!(ax, Float64[], Float64[], color=:red, label="A nullcline")
contour!(ax, xx, yy, ∂B416, levels=[0], color=:blue)
lines!(ax, Float64[], Float64[], color=:blue, label="B nullcline")
for sol in sols
    lines!(ax, sol, idxs=(1, 2), color=:tomato, label=nothing)
end
axislegend(ax, position = :rc)
limits!(ax, 0.0, 4.0, 0.0, 4.0)
fig

# ## Fig 4.17
prob417 = remake(prob415, p=ps416, u0=CA.ComponentArray(A=2.0, B=1.5), tspan=(0.0, 10.0))
@time sol = solve(prob417, Tsit5())

fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.17",
    aspect = 1,
)

lines!(ax, sol, idxs=(1, 2), color=:tomato, label="Trajectory")
xx = 1:0.01:3
yy = 1:0.01:3
∂A417 = [_dA415(x, y, ps416, nothing) for x in xx, y in yy]
∂B417 = [_dB415(x, y, ps416, nothing) for x in xx, y in yy];

streamplot!(ax, ∂F416, xx, yy)
contour!(ax, xx, yy, ∂A417, levels=[0], color=:red)
lines!(ax, Float64[], Float64[], color=:red, label="A nullcline")
contour!(ax, xx, yy, ∂B417, levels=[0], color=:blue)
lines!(ax, Float64[], Float64[], color=:blue, label="B nullcline")
axislegend(ax, position = :rc)
limits!(ax, 1.0, 3.0, 1.0, 3.0)
fig
