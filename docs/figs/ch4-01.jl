#===
# Chapter 4

Steady states and phase plots in an asymmetric network.
===#
using DifferentialEquations
using SimpleUnPack              ## @unpack
using Plots
Plots.default(linewidth=2)

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

# The model and derivatives for figure 4.1, 4.2, and 4.3.

function ∂A41(u, p, t)
    a, b = u
    @unpack k1, k2, k3, k4, k5, n = p
    return da = k1 * hil(1, b, n) - (k3 + k5) * a
end

function ∂B41(u, p, t)
    a, b = u
    @unpack k1, k2, k3, k4, k5, n = p
    return db = k2 + k5 * a - k4 * b
end

function model41!(D, u, p, t)
    D[1] = ∂A41(u, p, t)
    D[2] = ∂B41(u, p, t)
    return nothing
end

# ## Fig 4.2 A

ps1 = (k1=20., k2=5., k3=5., k4=5., k5=2., n=4.)
u0s = (
    [0.0, 0.0],
    [0.5, 0.6],
    [0.17, 1.1],
    [0.25, 1.9],
    [1.85, 1.7],
)

tend = 1.5

#---
sols = [solve(ODEProblem(model41!, u0, tend, ps1)) for u0 in u0s];

#---
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2 A (Time series)", labels=["[A]" "[B]"])

# ## Fig. 4.2 B (Phase plot)
plot( sols[1], idxs=(1, 2),
    xlabel="[A]", ylabel="[B]",
    aspect_ratio=:equal, legend=nothing,
    title="Fig. 4.2 B (Phase plot)",
    ylims=(0.0, 2.0), xlims=(0.0, 2.0)
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

plot!(fig, xlabel="[A]", ylabel="[B]", xlims=(0., 2.), ylims=(0., 2.), aspect_ratio=:equal)

# Let's sketch vector fields in phase plots.
∂A = function (x, y)
    ∂A41((x, y), ps1, 0.0)
end

∂B = function (x, y)
    ∂B41((x, y), ps1, 0.0)
end

∂F44 = function (x, y; scale=20)
    da = ∂A(x, y)
    db = ∂B(x, y)
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
    plot!(fig, sol, idxs=(1, 2), linealpha=0.7, lab=nothing)
end

## nullclines
contour!(fig, 0:0.01:2, 0:0.01:2, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")  ## Adding a fake line for the legend of A nullcline
contour!(fig, 0:0.01:2, 0:0.01:2, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline") ## Adding a fake line for the legend of B nullcline
plot!(fig, xlim=(0., 2.), ylim=(0., 2.), legend=:bottomright, size=(600, 600), xlabel="[A]", ylabel="[B]", aspect_ratio=:equal)

# ## Figure 4.5 B

# Vector field
fig = plot(title="Fig. 4.5 B (Vector field with nullclines)")
quiver!(fig, xx, yy, quiver=∂F44, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)

# Nullclines
contour!(fig, 0:0.01:2, 0:0.01:2, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")  ## Adding a fake line for the legend of A nullcline
contour!(fig, 0:0.01:2, 0:0.01:2, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline") ## Adding a fake line for the legend of B nullcline
plot!(fig, xlim=(0., 2.), ylim=(0., 2.), legend=:bottomright, size=(600, 600), xlabel="[A]", ylabel="[B]")
