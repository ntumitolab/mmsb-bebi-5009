#===
# Fig 6.07

Goldbeter Koshland switch
===#
using Plots
Plots.default(linewidth=2)

f607(w, K1, K2) = w * (1 - w + K1)/((1 - w) * (w + K2))

fig = plot(title = "Fig 6.07", xlims=(0, 3), ylims=(0, 1))
plot!(fig, x -> f607(x, 1, 1), identity, 0, 1, label="K1=K2=1")
plot!(fig, x -> f607(x, 0.1, 0.1), identity, 0, 1, label="K1=K2=0.1")
plot!(fig, x -> f607(x, 0.01, 0.01), identity, 0, 1, label="K1=K2=0.01")
