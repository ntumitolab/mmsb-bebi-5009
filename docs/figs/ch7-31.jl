# # Fig. 7.31
# Hasty synthetic oscillator model
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

#---
function model731!(du, u, p, t)
    @unpack alpha, sigma, gammax, gammay, ay, D = p
    @unpack x1, y1, x2, y2 = u
    den1 = (1 + x1^2 + sigma * x1^4) * (1 + y1^4)
    den2 = (1 + x2^2 + sigma * x2^4) * (1 + y2^4)
    du.x1 = (1 + x1^2 + alpha * sigma * x1^4) / den1 - gammax * x1 + D * (x2 - x1)
    du.y1 = ay * (1 + x1^2 + alpha * sigma * x1^4) / den1 - gammay * y1
    du.x2 = (1 + x2^2 + alpha * sigma * x2^4) / den2 - gammax * x2 + D * (x1 - x2)
    du.y2 = ay * (1 + x2^2 + alpha * sigma * x2^4) / den2 - gammay * y2
    nothing
end

ps731 = ComponentArray(
    alpha=11,
    sigma=2,
    gammax=0.2,
    gammay=0.012,
    ay=0.2,
    D=0.015,
)

u0731 = ComponentArray(
    x1=0.3963,
    y1=2.3346,
    x2=0.5578,
    y2=1.9317,
)

tend = 500.0
prob731 = ODEProblem(model731!, u0731, (0.0, tend), ps731)
@time sol731 = solve(prob731, KenCarp47())

#---
fig731 = Figure()
ax731 = Axis(fig731[1, 1];
    xlabel="Time (a.u.)",
    ylabel="Concentration (a.u.)",
    title="Fig. 7.31\nHasty Synthetic Oscillator Model",
)
lines!(ax731, 0..tend, t->sol731(t).x1, label="x1")
lines!(ax731, 0..tend, t->sol731(t).y1, label="y1")
lines!(ax731, 0..tend, t->sol731(t).x2, label="x2")
lines!(ax731, 0..tend, t->sol731(t).y2, label="y2")
axislegend(ax731; position=:rc)
fig731
