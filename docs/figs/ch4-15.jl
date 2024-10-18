# # Fig 4.15, 4.16, 4.17
# Oscillatory networks.
# ## Figure 4.15 (A)
using OrdinaryDiffEq
using Plots
Plots.default(linewidth=2)

#---
function dA415(u, p, t)
    a, b = u
    k0, k1, k2, n = p
    return dA = k0 - k1 * a * (1 + b^n)
end

function dB415(u, p, t)
    a, b = u
    k0, k1, k2, n = p
    return dB = k1 * a * (1 + b^n) - k2 * b
end

function model415!(D, u, p, t)
    D[1] = dA415(u, p, t)
    D[2] = dB415(u, p, t)
    return nothing
end

#---
ps1 = (k0 = 8., k1 = 1., k2 = 5., n = 2.)
u0s = (
    [1.5, 1.0],
    [0.0, 1.0],
    [0.0, 3.0],
    [2.0, 0.0],
)

tend = 8.
sols = map(u0s) do u0
    solve(ODEProblem(model415!, u0, tend, ps1))
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title ="Fig 4.15 (A)", xlims=(0., 8.))

# ## Fig 4.15 (B)

∂F415 = function (x, y; scale=20)
    da = dA415((x, y), ps1, 0.0)
    db = dB415((x, y), ps1, 0.0)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end

∂A = function (x, y)
    dA415((x, y), ps1, 0.0)
end

∂B = function (x, y)
    dB415((x, y), ps1, 0.0)
end

# Grid points
r = range(0, 4, 21)
xx = [x for y in r, x in r]
yy = [y for y in r, x in r];

#---
fig = plot(title="Fig 4.15 B")
quiver!(fig, xx, yy, quiver=∂F415, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 0:0.01:4, 0:0.01:4, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 0:0.01:4, 0:0.01:4, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
for sol in sols
    plot!(fig, sol, idxs=(1, 2), label=nothing)
end
plot!(fig, xlim=(0, 4), ylim=(0, 4), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

#===
## Fig 4.16 A

Oscillatory parameter set
===#

ps2 = (k0 = 8., k1 = 1., k2 = 5., n = 2.5)
tend = 100.0
u0s = (
    [1.5, 1.0],
    [0.0, 1.0],
    [0.0, 3.0],
    [2.0, 0.0],
)

sols = map(u0s) do u0
    solve(ODEProblem(model415!, u0, tend, ps2))
end;

plot(sols[1], xlabel="Time", ylabel="Concentration", title ="Fig 4.16(A)", xlims=(0., 8.))

# ## Fig 4.16 b
∂F416 = function (x, y; scale=20)
    da = dA415((x, y), ps2, 0.0)
    db = dB415((x, y), ps2, 0.0)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end

∂A = function (x, y)
    dA415((x, y), ps2, 0.0)
end

∂B = function (x, y)
    dB415((x, y), ps2, 0.0)
end

r = range(0, 4, 21)
xx = [x for y in r, x in r]
yy = [y for y in r, x in r];

#---
fig = plot(title="Fig 4.16 B")
quiver!(fig, xx, yy, quiver=∂F416, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 0:0.01:4, 0:0.01:4, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 0:0.01:4, 0:0.01:4, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
for sol in sols
    plot!(fig, sol, idxs=(1, 2), label=nothing)
end
plot!(fig, xlim=(0, 4), ylim=(0, 4), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

# ## Fig 4.17

sol = solve(ODEProblem(model415!, [2.0, 1.5], 10.0, ps2))

fig = plot(title="Fig 4.17")
plot!(fig, sol, idxs=(1, 2), label=nothing, arrow=:closed)
quiver!(fig, xx, yy, quiver=∂F416, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 1:0.01:3, 1:0.01:3, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, identity, 0, 0, line=(:black), label="A nullcline")
contour!(fig, 1:0.01:3, 1:0.01:3, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, identity, 0, 0, line=(:black, :dash), label="B nullcline")
plot!(fig, xlims=(1, 3), ylims=(1, 3), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")
