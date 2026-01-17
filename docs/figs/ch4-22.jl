#===
# Fig 4.22

Tangent line.
===#
using CairoMakie

# The curve
f = x -> 3 / (x-2)

#---
fig = Figure()
ax = Axis(fig[1, 1],
    title = "Fig 4.22",
    xlabel = "Reaction rate",
    ylabel = "Inhibitor concentration"
)
lines!(ax, 2.2..8.0, f, label="Curve")
lines!(ax, 2.7..5.3, x -> -3 / (4 - 2)^2 * (x - 4) + f(4), label="Tangent line")
axislegend(ax)
fig
