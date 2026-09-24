# # Fig 6.18
# Model of calcium-induced calcium release in hepatocytes
using Catalyst
using OrdinaryDiffEq
using Plots

@time "Build system" rn618 = @reaction_network begin
    @discretes I(t)
    (k1 * I, km1), R <--> RI
    (k2 * C, km2), RI <--> RIC
    (k3 * C, km3), RIC <--> RICC
    vr * (γ0 + γ1 * RIC) * (Cer - C), 0 --> C
    hill(C, p1, p2, 4), C => 0
end

#---
alg = FBDF()
up618 = Dict(
    :C => 0.0,
    :R => 1.0,
    :RI => 0.0,
    :RIC => 0.0,
    :RICC => 0.0,
    :k1 => 12.0,
    :km1 => 8.0,
    :k2 => 15.0,
    :km2 => 1.65,
    :k3 => 1.8,
    :km3 => 0.21,
    :vr => 0.185,
    :γ0 => 0.1,
    :γ1 => 20.5,
    :p1 => 8.5,
    :p2 => 0.065,
    :Cer => 8.37,
    :I => 0.0
)

# ## Fig 6.18 (A)
tspan = (0.0, 25.0)
up618a = merge(up618, Dict(:I => 2.0))
@time "Build problem" prob618a = ODEProblem(rn618, up618a, tspan; remove_conserved=true)
@time "Solve problem" sol618a = solve(prob618a, alg)
plot(sol618a, idxs=[:C, :RIC, :RICC], title="Fig 6.18 (A)", xlabel="Time", ylabel="Abundance", legend=:topright)

# ## Fig 6.18 (B)
@unpack I = rn618
events = [
    ModelingToolkitBase.SymbolicDiscreteCallback([20.0] => [I ~ 0.7]; discrete_parameters=[I]),
    ModelingToolkitBase.SymbolicDiscreteCallback([60.0] => [I ~ 1.2]; discrete_parameters=[I]),
    ModelingToolkitBase.SymbolicDiscreteCallback([90.0] => [I ~ 4.0]; discrete_parameters=[I])
]

osys618 = ode_model(rn618; remove_conserved=true, discrete_events=events)
tend = 120.0
prob618b = ODEProblem(osys618 |> complete, up618, tend)
@time sol618b = solve(prob618b, alg)
plot(sol618b, idxs=[:C], title="Fig 6.18 (B)", label="Calcium", xlabel="Time", ylabel="Concentration", legend=:topright)
