#===
# Fig 3.13

Generalized mass action (GMA) vs. Michaelis-Menten rate laws
===#
using Plots
Plots.gr(framestyle = :box)

plot([t -> 2t / (1+t), t -> t^0.4], 0, 4, label=["MM" "GMA"], xlabel="Substrate concentration (AU)", ylabel="Reaction rate (AU)", title="Fig 3.13")
