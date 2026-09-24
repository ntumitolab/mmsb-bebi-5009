# # Fig 7.07
# Model of lac operon in E. coli
using Startup: meshgrid, get_gradient, normalize_gradient
using DifferentialEquations
using SteadyStateDiffEq
using DiffEqCallbacks
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

#---
@time "Build system" rn707 = @reaction_network begin
    (a1 / (1 + RToverK1 * (K2 / (K2 + L))^4), δM), 0 <--> M
    (c1 * M, δY), 0 <--> Y
    mm(Le, kL * Y, KML), 0 --> L
    δL, L --> 0
    (Y / 4) * mm(L, 2kg, KMg), L => 0
end

#---
up707 = Dict(
    :δM => 0.48, :δY => 0.03, :δL => 0.02,
    :a1 => 0.29, :K2 => 2.92e6, :RToverK1 => 213.2, :c1 => 18.8,
    :kL => 6e4, :KML => 680.0, :kg => 3.6e3, :KMg => 7e5,
    :Le => 0.0, :M => 0.01, :Y => 0.1, :L => 0.0,
)

@unpack Le = rn707
cbs = let
    event_le1 = PresetTimeCallback([500.0], (integrator) -> integrator.ps[Le] = 50)
    event_le2 = PresetTimeCallback([1000.0], (integrator) -> integrator.ps[Le] = 100)
    event_le3 = PresetTimeCallback([1500.0], (integrator) -> integrator.ps[Le] = 150)
    event_le4 = PresetTimeCallback([2000.0], (integrator) -> integrator.ps[Le] = 0)
    CallbackSet(event_le1, event_le2, event_le3, event_le4)
end

alg = FBDF()

# ## Fig 7.07 (A)
tend = 2500.0
@time "Build problem" prob707a = ODEProblem(rn707, up707, tend)
@time "Solve problem" sol = solve(prob707a, alg, callback=cbs)

plot(sol, idxs=[:Y], labels="β-galactosidase monomer", title="Fig 7.07 (A)", xlabel="Time (min)", ylabel="Concentration")
plot!(t -> 50 * (500 <= t < 1000) + 100 * (1000 <= t < 1500) + 150 * (1500 <= t < 2000), 0, tend, linestyle=:dash, label="External lactose")

# ## Fig 7.07 (B)
# Comparing the original model and the modified model
@time "Build system" rn707b = @reaction_network begin
    (a1 / (1 + RToverK1 * (K2 / (K2 + L))^4), δM), 0 <--> M
    (c1 * M, δY), 0 <--> Y
    Enz * mm(Le, 4kL, KML), 0 --> L
    δL, L --> 0
    Enz * mm(L, 2kg, KMg), L => 0
end

#---
up707b = merge(up707, Dict(:Enz => 40.0))
prob707a = SteadyStateProblem(rn707, up707)
prob707b = SteadyStateProblem(rn707b, up707b)
lerange = 0:100
prob_func = (prob, ctx) -> remake(prob, p=[:Le => lerange[ctx.sim_id]])
output_func = (sol, ctx) -> (sol[:Y] / 4, false)
ssalg = DynamicSS(FBDF())

eprob = EnsembleProblem(prob707a; prob_func, output_func)
eprob_mod = EnsembleProblem(prob707b; prob_func, output_func)
@time sim = solve(eprob, ssalg; trajectories=length(lerange))
@time sim_mod = solve(eprob_mod, ssalg; trajectories=length(lerange))

plot(lerange, [sim.u sim_mod.u], labels=["Original model" "Modified model"], title="Fig 7.07 (B)", xlabel="External lactose concentration", ylabel="β-galactosidase tetramer concentration")
