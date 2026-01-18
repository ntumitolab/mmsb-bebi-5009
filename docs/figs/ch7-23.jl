# # Fig 7.23
# Hasty synthetic oscillator model
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

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
@time sol723 = solve(prob723, KenCarp47())

#---
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig 7.23 A"
)
lines!(ax, 0..tend, t-> sol723(t).x, label = "x")
lines!(ax, 0..tend, t-> sol723(t).y, label = "y")
axislegend(ax, position = :rc)
fig

#---
xx = range(0, 1.5, 51)
yy = range(0, 3, 51)
tt = 0:1:tend
∂X723 = [_dx723((; x=x, y=y), ps723, nothing) for x in xx, y in yy]
∂Y723 = [_dy723((; x=x, y=y), ps723, nothing) for x in xx, y in yy]
aa = [sol723(t).x for t in tt]
bb = [sol723(t).y for t in tt];

# Vector field
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "X",
    ylabel = "Y",
    title = "Fig 7.23 B\nVector field",
    aspect = 1,
)
contour!(ax, xx, yy, ∂X723, levels=[0], color=:black, linestyle=:solid, linewidth=2, label="X nullcline")
contour!(ax, xx, yy, ∂Y723, levels=[0], color=:black, linestyle=:dash, linewidth=2, label="Y nullcline")
lines!(ax, aa, bb, color=:tomato, label="Trajectory")
axislegend(ax, position = :rb)
limits!(ax, 0.0, 1.5, 1.5, 3.0)
fig
