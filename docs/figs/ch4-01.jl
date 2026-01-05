#===
# Fig 4.2, 4.3, 4.4, 4.5

Steady states and phase plots in an asymmetric network.
===#
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# The model for figure 4.1, 4.2, and 4.3.
_dA401(u, p, t) = p.k1 / (1 + u.B^p.n) - (p.k3 + p.k5) * u.A
_dB401(u, p, t) = p.k2 + p.k5 * u.A - p.k4 * u.B

function model401!(D, u, p, t)
    D.A = _dA401(u, p, t)
    D.B = _dB401(u, p, t)
    return nothing
end

# ## Fig 4.2 A
tend = 1.5
ps402a = ComponentArray(
    k1 = 20.0,
    k2 = 5.0,
    k3 = 5.0,
    k4 = 5.0,
    k5 = 2.0,
    n = 4.0
)
ics402a = ComponentArray(
    A = 0.0,
    B = 0.0
)
prob401 = ODEProblem(model401!, ics402a, (0.0, tend), ps402a)

u0s = [
    ComponentArray(
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
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2 A (Time series)", labels=["A" "B"], legend=:right)

# ## Fig. 4.2 B (Phase plot)
plot(
    sols[1], idxs=(1, 2),
    xlabel="[A]", ylabel="[B]",
    aspect_ratio=:equal, legend=nothing,
    title="Fig. 4.2 B (Phase plot)",
    ylims=(0.0, 2.0), xlims=(0.0, 2.0),
    size=(600, 600)
)

# ## Fig. 4.3 A (Multiple time series)
fig = plot(title="Fig. 4.3A (Multiple time series)")

for sol in sols
    plot!(fig, sol, linealpha=0.5, legend=nothing)
end

plot!(fig, xlabel="Time", ylabel="Concentration")

# ## Fig. 4.3 B (Phase plot)
fig = plot(title="Fig. 4.3 B (Phase plot)")

for sol in sols
    plot!(fig, sol, idxs=(1, 2), legend=nothing)
end

plot!(fig, xlabel="[A]", ylabel="[B]", xlims=(0.0, 2.0), ylims=(0.0, 2.0), aspect_ratio=:equal, size=(600, 600))

# Let's sketch vector fields in phase plots.
∂F44 = function (x, y; scale=20)
    da = _dA401(ComponentArray(A=x, B=y), ps402a, nothing)
    db = _dB401(ComponentArray(A=x, B=y), ps402a, nothing)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end

# Grid points
rxy = range(0.0, 2.0, 11)
xx = [x for y in rxy, x in rxy]
yy = [y for y in rxy, x in rxy];

# ## Figure 4.4A
fig = plot(title="Fig. 4.4 A (Phase plot with vector field)")

quiver!(fig, xx, yy, quiver=∂F44, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)

for sol in sols
    plot!(fig, sol, idxs=(1, 2), linealpha=0.7, legend=nothing)
end

plot!(fig, size=(600, 600), xlims=(rxy[1], rxy[end]), ylims=(rxy[1], rxy[end]), xlabel="[A]", ylabel="[B]")

# ## Figure 4.5A
fig = plot(title="Fig. 4.5 A (Phase plot with nullclines)")

## Phase plots
for sol in sols
    plot!(fig, sol, idxs=(1, 2), linealpha=0.7, label=false)
end

## nullclines
∂A44 = (x, y) -> _dA401(ComponentArray(A=x, B=y), ps402a, nothing)
∂B44 = (x, y) -> _dB401(ComponentArray(A=x, B=y), ps402a, nothing)
contour!(fig, 0:0.01:2, 0:0.01:2, ∂A44, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")  ## Adding a fake line for the legend of A nullcline
contour!(fig, 0:0.01:2, 0:0.01:2, ∂B44, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline") ## Adding a fake line for the legend of B nullcline
plot!(fig, xlim=(0.0, 2.0), ylim=(0.0, 2.0), legend=:bottomright, size=(600, 600), xlabel="[A]", ylabel="[B]", aspect_ratio=:equal)

# ## Figure 4.5 B
# Vector field
fig = plot(title="Fig. 4.5 B (Vector field with nullclines)")
quiver!(fig, xx, yy, quiver=∂F44, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)

# Nullclines
contour!(fig, 0:0.01:2, 0:0.01:2, ∂A44, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")  ## Adding a fake line for the legend of A nullcline
contour!(fig, 0:0.01:2, 0:0.01:2, ∂B44, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline") ## Adding a fake line for the legend of B nullcline
plot!(fig, xlim=(0.0, 2.0), ylim=(0.0, 2.0), legend=:bottomright, size=(600, 600), xlabel="[A]", ylabel="[B]")
