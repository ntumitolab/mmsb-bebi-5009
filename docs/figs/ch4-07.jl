# # Fig 4.7, 4.8
# Symmetric (bistable) biological networks.

using DifferentialEquations
using LabelledArrays
using UnPack
using Plots
Plots.default(linewidth=2)
import DisplayAs.SVG

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

#---
function ∂A47(u, p, t)
    a, b = u
    @unpack k1, k2, k3, k4, n1, n2 = p
    return da = k1 * hil(1, b, n1) - k3 * a
end

function ∂B47(u, p, t)
    a, b = u
    @unpack k1, k2, k3, k4, n1, n2 = p
    return db = k2 * hil(1, a, n2) - k4 * b
end

function model47!(D, u, p, t)
    D.a = ∂A47(u, p, t)
    D.b = ∂B47(u, p, t)
    return nothing
end

#===
## Fig 4.7 A

Asymmetric parameter set
===#

ps1 = (k1=20., k2=20., k3=5., k4=5., n1=4., n2=1.)
tend = 4.0

sol1 = solve(ODEProblem(model47!, LVector(a=3., b=1.), tend, ps1))
sol2 = solve(ODEProblem(model47!, LVector(a=1., b=3.), tend, ps1))

ax1 = plot(sol1, xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.7A (1)")
ax2 = plot(sol2, xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.7A (2)")
fig = plot(ax1, ax2, layout=(2, 1), size=(600, 600))

fig |> SVG

# ## Fig 4.7 B

∂F47 = function (x, y; scale=20)
    da = ∂A47((x, y), ps1, 0.0)
    db = ∂B47((x, y), ps1, 0.0)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end

∂A = function (x, y)
    ∂A47((x, y), ps1, 0.0)
end

∂B = function (x, y)
    ∂B47((x, y), ps1, 0.0)
end

# Grid points
r = range(0., 5., 21)
xx = [x for y in r, x in r]
yy = [y for y in r, x in r];

fig = plot(title="Fig 4.7 B (Vector field with nullclines)")
quiver!(fig, xx, yy, quiver=∂F47, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)

contour!(fig, 0:0.01:5, 0:0.01:5, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 0:0.01:5, 0:0.01:5, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
plot!(fig, xlim=(0, 5), ylim=(0, 5), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

fig |> SVG

#===
## Fig 4.8

Symmetric parameter set
===#

ps2 = LVector(k1=20., k2=20., k3=5., k4=5., n1=4., n2=4.)

tend = 4.0
sol1 = solve(ODEProblem(model47!, LVector(a=3., b=1.), tend, ps2))
sol2 = solve(ODEProblem(model47!, LVector(a=1., b=3.), tend, ps2))

ax1 = plot(sol1, xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.8A (1)")
ax2 = plot(sol2, xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.8A (2)")
fig = plot(ax1, ax2, layout=(2, 1), size=(600, 600))

fig |> SVG


#---
∂F48 = function (x, y; scale=20)
    da = ∂A47((x, y), ps2, 0.0)
    db = ∂B47((x, y), ps2, 0.0)
    s = sqrt(hypot(da, db)) * scale
    return (da / s, db / s)
end

∂A = function (x, y)
    ∂A47((x, y), ps2, 0.0)
end

∂B = function (x, y)
    ∂B47((x, y), ps2, 0.0)
end

#---
fig = plot(title="Fig 4.8 B")
quiver!(fig, xx, yy, quiver=∂F48, line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)

contour!(fig, 0:0.01:5, 0:0.01:5, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 0:0.01:5, 0:0.01:5, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
plot!(fig, xlim=(0, 5), ylim=(0, 5), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

fig |> SVG

#===
## Fig 4.8 C
around the unstable steady-state
===#

r2 = range(1.0, 1.5, 16)
xx2 = [x for y in r2, x in r2]
yy2 = [y for y in r2, x in r2]

fig = plot(title="Fig 4.8 C (close up)")

quiver!(fig, xx2, yy2, quiver=(x, y)-> ∂F48(x, y; scale=60), line=(:lightgrey), arrow=(:closed), aspect_ratio=:equal)
contour!(fig, 1:0.01:1.5, 1:0.01:1.5, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="A nullcline")
contour!(fig, 1:0.01:1.5, 1:0.01:1.5, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="B nullcline")
plot!(fig, xlim=(1, 1.5), ylim=(1, 1.5), legend=:topright, size=(600, 600), xlabel="[A]", ylabel="[B]")

fig |> SVG

# Another way to draw nullclines is to find the analytical solution when dA (or dB) is zero. And then sketch the nullclines in a parameteric plot.

nca47(b, p) = p.k1 / p.k3 * hil(1, b, p.n1)
ncb47(a, p) = p.k2 / p.k4 * hil(1, a, p.n2)

pls = map((8.0, 16.0, 20.0, 35.0)) do k1
    ps = LVector(k1=k1, k2=20., k3=5., k4=5., n1=4., n2=4.)
    pl = plot(b -> nca47(b, ps), identity, 0., 7., label="Nullcline A")
    plot!(pl, identity, a -> ncb47(a, ps), 0., 7., label="Nullcline B")
    plot!(pl, title = "K1 = $k1", xlim=(0., 7.), ylim=(0., 7.), aspect_ratio = :equal, xlabel="[A]", ylabel="[B]")
    pl
end

fig = plot(pls..., size = (800, 800))

fig |> SVG

# ## Runtime information
import Pkg
Pkg.status()

#---
import InteractiveUtils
InteractiveUtils.versioninfo()
