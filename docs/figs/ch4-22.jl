#===
# Fig 4.22

Tangent line.

Using ForwardDiff to calculate the derivative on a curve.
===#

using ForwardDiff
using Plots
Plots.default(linewidth=2)
import DisplayAs.SVG

# The curve
f = x -> 3 / (x-2)

# slope and y intercept for x = 4
slope = ForwardDiff.derivative(f, 4)
yintercept = f(4)

#---
g = x -> slope * (x - 4) + yintercept

fig = plot(title="Fig 4.22")
plot!(fig, f, 2.2, 8.0, lab="Curve")
plot!(fig, g, 2.7, 5.3, lab="Tangent line")
plot!(fig, xlabel="Reaction rate", ylabel="Inhibitor concentration",
      xlims=(2., 8.), ylims=(0., 4.)
)

fig |> SVG
