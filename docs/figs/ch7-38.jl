# # Fig. 7.38
# stochastic implementation (Gillespie SSA) of a collection of spontaneously decaying molecules.
using JumpProcesses
using Catalyst
using Plots

@time "Build system" rn738 = @reaction_network begin
    d, A --> 0
end

tend = 5.0
ps738 = Dict(:d => 1.0)
@time "Build problem" jprob1000 = JumpProblem(rn738, [:A => 1000], (0.0, tend), ps738)
jprob100 = remake(jprob1000, u0=[:A => 100])
jprob10 = remake(jprob1000, u0=[:A => 10])

@time "Solve problem" sol1000 = solve(jprob1000)
@time "Solve problem" sol100 = solve(jprob100)
@time "Solve problem" sol10 = solve(jprob10)

#---
plot(sol1000, xlabel="Time", ylabel="# of molecules", title="Fig. 7.38 (A)", label="1000 molecules")
plot!(t -> 1000 * exp(-t), linestyle=:dash, label="ODE solution")

#---
plot(sol100, xlabel="Time", ylabel="# of molecules", title="Fig. 7.38 (B)", label="100 molecules")
plot!(t -> 100 * exp(-t), linestyle=:dash, label="ODE solution")

#---
plot(sol10, xlabel="Time", ylabel="# of molecules", title="Fig. 7.38 (C)", label="10 molecules")
plot!(t -> 10 * exp(-t), linestyle=:dash, label="ODE solution")
