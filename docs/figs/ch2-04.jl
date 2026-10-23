# Fig 2.04
# Exponential decay
using Plots
Plots.default(linewidth=1.5)

pl = plot(xlabel="Time", ylabel="Concentration", title="Fig 2.4")
for k in 1:3
    plot!(pl, t -> 3 * exp(-k * t), 0, 5, label="exp(-$(k)t)")
end
pl
