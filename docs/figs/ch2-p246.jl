# # Prob 2.4.6
using OrdinaryDiffEq
using CairoMakie

#---
sol = ODEProblem((u, p, t) -> p * (1.0 - u), 0.0, 10.0, 1.0) |> solve

# Fig 2.46
fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="Concentration",
    title="Prob. 2.4.6"
)
lines!(ax, 0 .. 10.0, t -> sol(t), label="u(t)")
axislegend(ax, position=:rb)
fig
