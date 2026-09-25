# # Fig 4.22
# Tangent line.
using Plots
using ForwardDiff
Plots.gr(linewidth=1.5)

f = x -> 3 / (x-2)
g = x -> ForwardDiff.derivative(f, 4) * (x - 4) + f(4)
plot(f, 2.2, 8.0, label="f(x)")
plot!(g, 2.7, 5.3, label="Tangent line")
plot!(title="Fig. 4.22", xlabel="x", ylabel="f(x)", legend=:topright, xlims=(2, 6), ylims=(0, 5))
