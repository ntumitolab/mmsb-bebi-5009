# # Fig 4.15, 4.16, 4.17
# Oscillatory networks.
# ## Figure 4.15 (A)
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# The model
_dA415(u, p, t) = p.k0 - p.k1 * u.A * (1 + u.B^p.n)
_dB415(u, p, t) = p.k1 * u.A * (1 + u.B^p.n) - p.k2 * u.B
function model415!(D, u, p, t)
    D.A = _dA415(u, p, t)
    D.B = _dB415(u, p, t)
    nothing
end

#---
ps415 = ComponentArray(
    k0 = 8.0,
    k1 = 1.0,
    k2 = 5.0,
    n = 2.0
)
u0415 = ComponentArray(
    A = 1.5,
    B = 1.0
)

tend = 8.0
prob415 = ODEProblem(model415!, u0415, (0.0, tend), ps415)

#---
u0s = [
    ComponentArray(A=1.5, B=1.0),
    ComponentArray(A=0.0, B=1.0),
    ComponentArray(A=0.0, B=3.0),
    ComponentArray(A=2.0, B=0.0),
]

@time sols = map(u0s) do u0
    solve(remake(prob415, u0=u0))
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title ="Fig 4.15 (A)", xlims=(0, 8), labels=["A" "B"])

# ## Fig 4.15 (B)

∂F415 = function (x, y; scale=20)
    da = _dA415(ComponentArray(A=x, B=y), ps415, nothing)
    db = _dB415(ComponentArray(A=x, B=y), ps415, nothing)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end

∂A415 = (x, y) -> _dA415(ComponentArray(A=x, B=y), ps415, nothing)
∂B415 = (x, y) -> _dB415(ComponentArray(A=x, B=y), ps415, nothing)

# Grid points
r = range(0, 4, 21)
xx = [x for y in r, x in r]
yy = [y for y in r, x in r];

#---
fig = plot(title="Fig 4.15 B")
quiver!(fig, xx, yy, quiver=∂F415, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 0:0.01:4, 0:0.01:4, ∂A415, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 0:0.01:4, 0:0.01:4, ∂B415, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
for sol in sols
    plot!(fig, sol, idxs=(1, 2), label=nothing)
end
plot!(fig, xlim=(0, 4), ylim=(0, 4), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

#===
## Fig 4.16 A

Oscillatory parameter set
===#
ps416 = ComponentArray(
    k0 = 8.0,
    k1 = 1.0,
    k2 = 5.0,
    n = 2.5
)

tend = 100.0
u0s = [
    ComponentArray(A=1.5, B=1.0),
    ComponentArray(A=0.0, B=1.0),
    ComponentArray(A=0.0, B=3.0),
    ComponentArray(A=2.0, B=0.0),
]

prob416 = remake(prob415, p=ps416)

#---
@time sols = map(u0s) do u0
    solve(remake(prob416, u0=u0))
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title ="Fig 4.16(A)", xlims=(0, 8), labels=["A" "B"])

# ## Fig 4.16 b
∂F416 = function (x, y; scale=20)
    da = _dA415(ComponentArray(A=x, B=y), ps416, nothing)
    db = _dB415(ComponentArray(A=x, B=y), ps416, nothing)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end
∂A416 = (x, y) -> _dA415(ComponentArray(A=x, B=y), ps416, nothing)
∂B416 = (x, y) -> _dB415(ComponentArray(A=x, B=y), ps416, nothing)

#---
r = range(0, 4, 21)
xx = [x for y in r, x in r]
yy = [y for y in r, x in r];

#---
fig = plot(title="Fig 4.16 B")
quiver!(fig, xx, yy, quiver=∂F416, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 0:0.01:4, 0:0.01:4, ∂A416, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 0:0.01:4, 0:0.01:4, ∂B416, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
for sol in sols
    plot!(fig, sol, idxs=(1, 2), label=nothing)
end
plot!(fig, xlim=(0, 4), ylim=(0, 4), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

# ## Fig 4.17
prob417 = remake(prob415, p=ps416, u0=ComponentArray(A=2.0, B=1.5), tspan=(0.0, 10.0))
sol = solve(prob417)

fig = plot(title="Fig 4.17")
plot!(fig, sol, idxs=(1, 2), label=nothing, arrow=:closed)
quiver!(fig, xx, yy, quiver=∂F416, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 1:0.01:3, 1:0.01:3, ∂A416, levels=[0], cbar=false, line=(:black))
plot!(fig, identity, 0, 0, line=(:black), label="A nullcline")
contour!(fig, 1:0.01:3, 1:0.01:3, ∂B416, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, identity, 0, 0, line=(:black, :dash), label="B nullcline")
plot!(fig, xlims=(1, 3), ylims=(1, 3), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")
