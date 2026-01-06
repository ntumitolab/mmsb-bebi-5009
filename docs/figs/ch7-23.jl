# # Fig 7.23
# Hasty synthetic oscillator model

using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

#---
_dx723(u, p, t) = (1 + u.x^2 + p.alpha * p.sigma * u.x^4) / ((1 + u.x^2 + p.sigma * u.x^4) * (1 + u.y^4)) - p.gammax * u.x
_dy723(u, p, t) = p.ax * ((1 + u.x^2 + p.alpha * p.sigma * u.x^4) / ((1 + u.x^2 + p.sigma * u.x^4) * (1 + u.y^4))) - p.gammay * u.y
function model723!(D, u, p, t)
    @unpack alpha, sigma, gammax, gammay, ax = p
    @unpack x, y = u
    D.x = _dx723(u, p, t)
    D.y = _dy723(u, p, t)
    nothing
end

#---
ps723 = ComponentArray(
    alpha=11.0,
    sigma=2.0,
    gammax=0.2,
    gammay=0.012,
    ax=0.2,
)

ics723 = ComponentArray(
    x=0.3963,
    y=2.3346,
)

tend = 300.0
prob723 = ODEProblem(model723!, ics723, (0.0, tend), ps723)

#---
@time sol723 = solve(prob723, Tsit5())

#---
plot(sol723, xlabel="Time", ylabel="Concentration", title="Fig 7.23 A", labels=["x" "y"])

#---
∂X723 = (x, y) -> _dx723(ComponentArray(x=x, y=y), ps723, nothing)
∂Y723 = (x, y) -> _dy723(ComponentArray(x=x, y=y), ps723, nothing)

# Grid points
rX = range(0, 1.5, 31)
rY = range(0, 3, 31)
xx = [x for y in rY, x in rX]
yy = [y for y in rY, x in rX];

# Vector field
fig = plot(title="Fig 7.23 B")
contour!(fig, rX, rY, ∂X723, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="X nullcline")
contour!(fig, rX, rY, ∂Y723, levels=[0], cbar=false, line=(:dot, :black))
plot!(fig, Float64[], Float64[], line=(:dot, :black), label="Y nullcline")
plot!(sol723, idxs=(1,2), lw=2, label=false)
plot!(fig, xlims=(0.0, 1.5), ylims=(1.5, 3.0), xlabel="X", ylabel="Y")
