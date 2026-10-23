# Fig 3.03
# Michaelis-Menten kinetics
using OrdinaryDiffEq
using Catalyst
using Plots

# Enzyme kinetics full model
@time "Build system" rn303 = @reaction_network begin
    k1, S + E --> ES
    km1, ES --> S + E
    k2, ES --> P + E
end

#---
up303 = Dict(:S => 5.0, :ES => 0.0, :P => 0.0, :E=>1.0, :k1 => 30.0, :km1 => 1.0, :k2 => 10.0)
tend = 1.0
@time "Build problem" prob303 = ODEProblem(rn303, up303, (0.0, tend); remove_conserved=true)
@time "Solve problem" sol303 = solve(prob303, Tsit5())
plot(sol303, idxs=[:S, :ES, :P, :E], title="Fig. 3.03 (Full model)", xlabel="Time", ylabel="Concentration")

# ## QSSA of the ES complex
@time "Build system" rn303mm = @reaction_network begin
    mm(S, k2 * ET, (km1 + k2) / k1), S => P
end

#---
up303mm = Dict(:k1 => 30.0, :km1 => 1.0, :k2 => 10.0, :ET => 1.0, :S => 5.0, :P => 0.0)
tend = 1.0
@time "Build problem" prob303mm = ODEProblem(rn303mm, up303mm, (0.0, tend); remove_conserved=true)
@time "Solve problem" sol303mm = solve(prob303mm, Tsit5())
plot(sol303mm, idxs=[:S, :P], labels=["S (QSSA)" "P (QSSA)"], title="Fig. 3.03 (QSSA vs Full)", xlabel="Time", ylabel="Concentration")
plot!(sol303, idxs=[:S, :P], labels=["S (full)" "P (full)"], linestyle=:dash, legend=:right)
