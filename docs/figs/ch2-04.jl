#===
# Fig 2.04
Exponential decay
===#

using Plots
Plots.default(linewidth=2)
import DisplayAs.SVG

fig = plot(title= "Fig 2.4")
for k in 1:3
    plot!(fig, t -> 3 * exp(-k*t), 0., 5., label = "exp(-$(k)t)")
end

plot!(fig, xlim = (0, 5), ylim=(0, 3.2), xlabel="Time", ylabel="Concentration")

fig |> SVG
