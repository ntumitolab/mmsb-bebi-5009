# # Fig. 7.38
# stochastic implementation (Gillespie SSA) of a collection of spontaneously decaying molecules.
using CairoMakie
using JumpProcesses
using ComponentArrays: ComponentArray

# Define propensity functions and effects
v1_738(u, p, t) = p.d * u[1]
function affectv1_738!(integrator)
    integrator.u[1] -= 1        ## A -> 01
    nothing
end
#---
ps738 = ComponentArray(d=1)
tend = 5.0
jump1 = ConstantRateJump(v1_738, affectv1_738!)

prob1000 = DiscreteProblem(ComponentArray(A=1000), (0.0, tend), ps738)
jump_prob1000 = JumpProblem(prob1000, Direct(), jump1)
@time sol1000 = solve(jump_prob1000, SSAStepper())

prob100 = DiscreteProblem(ComponentArray(A=100), (0.0, tend), ps738)
jump_prob100 = JumpProblem(prob100, Direct(), jump1)
@time sol100 = solve(jump_prob100, SSAStepper())

prob10 = DiscreteProblem(ComponentArray(A=10), (0.0, tend), ps738)
jump_prob10 = JumpProblem(prob10, Direct(), jump1)
@time sol10 = solve(jump_prob10, SSAStepper())

#---
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel="Time",
    ylabel="# of molecules",
    title="A"
)
lines!(ax, sol1000.t, reduce(vcat, sol1000.u))
lines!(ax, 0..tend, t->1000 * exp(-ps738.d * t), linestyle=:dash)
limits!(ax, 0, tend, 0, 1000)

ax = Axis(fig[1, 2],
    xlabel="Time",
    ylabel="# of molecules",
    title="B"
)
lines!(ax, sol100.t, reduce(vcat, sol100.u))
lines!(ax, 0..tend, t->100 * exp(-ps738.d * t), linestyle=:dash)
limits!(ax, 0, tend, 0, 100)

ax = Axis(fig[1, 3],
    xlabel="Time",
    ylabel="# of molecules",
    title="C"
)
lines!(ax, sol10.t, reduce(vcat, sol10.u))
lines!(ax, 0..tend, t->10 * exp(-ps738.d * t), linestyle=:dash)
limits!(ax, 0, tend, 0, 10)
fig
