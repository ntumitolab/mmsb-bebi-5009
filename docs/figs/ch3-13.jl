#===
# Fig 3.13

Generalized mass action (GMA) vs. Michaelis-Menten rate laws
===#

using Plots
Plots.default(linewidth=2)

plot(
    [t -> 2t / (1+t), t -> t^0.4], 0., 4.,
    label = ["MM" "GMA"], title = "Fig 3.13",
    xlabel= "Substrate concentration (AU)",
    ylabel="Reaction rate (AU)"
)

# ## Runtime information
import Pkg
Pkg.status()

#---
import InteractiveUtils
InteractiveUtils.versioninfo()
