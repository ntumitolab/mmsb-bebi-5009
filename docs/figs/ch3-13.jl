#===
# Fig 3.13

Generalized mass action (GMA) vs. Michaelis-Menten rate laws
===#
using CairoMakie

mm = t -> 2t / (1+t)
gma = t -> t^0.4

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Substrate concentration (AU)",
    ylabel = "Reaction rate (AU)",
    title = "Fig 3.13",
)
lines!(ax, 0..4, mm, label = "MM")
lines!(ax, 0..4, gma, label = "GMA")
axislegend(ax, position = :rb)
fig
