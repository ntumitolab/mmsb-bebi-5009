# # Fig 3.13
# GMA and Michaelis-Menten rate laws

import DisplayAs.SVG
using Plots
Plots.default(linewidth=2)

fig = plot(
    [t -> 2t / (1+t), t -> t^0.4], 0., 4.,
    label = ["MM" "GMA"], title = "Fig 3.13",
    xlabel= "Substrate concentration (AU)",
    ylabel="Reaction rate (AU)"
)

fig |> SVG

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()
