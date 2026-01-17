#===
# Fig 6.07

Goldbeter Koshland switch
===#
using CairoMakie

f607(w, K1, K2) = w * (1 - w + K1)/((1 - w) * (w + K2))
yy = 0:0.001:1
xx1 = [f607(y, 1, 1) for y in yy]
xx2 = [f607(y, 0.1, 0.1) for y in yy]
xx3 = [f607(y, 0.01, 0.01) for y in yy]

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Stimulus", ylabel="Response", title="Fig 6.07")
lines!(ax, xx1, yy, label="K1=K2=1")
lines!(ax, xx2, yy, label="K1=K2=0.1")
lines!(ax, xx3, yy, label="K1=K2=0.01")
limits!(ax, 0, 3, 0, 1)
axislegend(ax, position=:rc)
fig
