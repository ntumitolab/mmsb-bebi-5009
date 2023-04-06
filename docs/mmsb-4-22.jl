# # Figure 4.22
# Tangent line
import DisplayAs.PNG
using Plots
Plots.default(linewidth=2)

fig = plot(title="Fig 4.22 Tangent line")
plot!(fig, t -> 3 / (t-2), 2.2, 8.0, lab="Curve")
plot!(fig, t -> 1.5 - (t - 4) * 0.75, 2.7, 5.3, lab="Tangent line")
plot!(fig, xlabel="Reaction rate", ylabel="Inhibitor concentration",
      xlims=(2., 8.), ylims=(0., 4.))

fig |> PNG

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()
