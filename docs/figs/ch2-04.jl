#===
# Fig 2.04
Exponential decay
===#
using CairoMakie

#---
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Fig 2.4\nExponential Decay"
)

for k in 1:3
    lines!(ax, 0 .. 5, t -> 3 * exp(-k * t), label="exp(-$(k)t)")
end

axislegend(ax, position=:rt)

fig
