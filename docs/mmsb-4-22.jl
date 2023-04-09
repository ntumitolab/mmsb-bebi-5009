# # Fig 4.22
# Tangent line
using ForwardDiff
using Plots
Plots.default(linewidth=2)

# Curve
f = x -> 3 / (x-2)

# slopeand y intercept for x = 4
slope = ForwardDiff.derivative(f, 4)
yintercept = f(4)
# Tangent line
g = x -> slope * (x - 4) + yintercept

fig = plot(title="Fig 4.22 Tangent line")
plot!(fig, f, 2.2, 8.0, lab="Curve")
plot!(fig, g, 2.7, 5.3, lab="Tangent line")
plot!(fig, xlabel="Reaction rate", ylabel="Inhibitor concentration",
      xlims=(2., 8.), ylims=(0., 4.))

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()
