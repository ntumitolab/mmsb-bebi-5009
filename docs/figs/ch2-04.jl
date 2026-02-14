#===
# Fig 2.04
Exponential decay
===#
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)
#---
pl1 = plot(xlabel="Time", ylabel="Concentration", title="Fig 2.4\nExponential Decay")

for k in 1:3
    plot!(pl1, t -> 3 * exp(-k * t), 0, 5, label="exp(-$(k)t)")
end

pl1
