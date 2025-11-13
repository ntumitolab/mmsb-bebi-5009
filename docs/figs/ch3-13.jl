#===
# Fig 3.13

Generalized mass action (GMA) vs. Michaelis-Menten rate laws
===#
using Plots
Plots.default(linewidth=2)

mm = t -> 2t / (1+t)
gma = t -> t^0.4

fig = plot(
    [mm, gma], 0., 4.,
    label = ["MM" "GMA"], title = "Fig 3.13",
    xlabel= "Substrate concentration (AU)",
    ylabel="Reaction rate (AU)"
)
