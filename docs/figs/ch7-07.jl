#===
# Fig 7.7

model of lac operon in E. coli
===#

using Catalyst
using ModelingToolkit
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

# PNG output in Literate.jl
PNG(fig) = display("image/png", fig)

#---
rn = @reaction_network begin
    a1 / (1 + RToverK1 * (K2 / (K2 + L))^4), 0 --> M
    δM, M --> 0
    (c1 * M, δY), 0 <--> Y
    mm(Le, kL * Y, KML), 0 --> L
    mm(L, 2 * kg * (Y / 4), KMg), L => 0
    δL, L --> 0
end

#---

setdefaults!(rn, [
    :δM => 0.48,
    :δY => 0.03,
    :δL => 0.02,
    :a1 => 0.29,
    :K2 => 2.92 * 1e6,
    :RToverK1 => 213.2,
    :c1 => 18.8,
    :kL => 6 * 1e4,
    :KML => 680,
    :kg => 3.6 * 1e3,
    :KMg => 7 * 1e5,
    :Le => 0.0,
    :M => 0.01,
    :Y => 0.1,
    :L => 0.0
])

osys = convert(ODESystem, rn; remove_conserved=true)

equations(osys)

# ## Fig 7.07 (A)

@unpack Le = osys
idx = findfirst(isequal(Le), parameters(osys))

cb1 = PresetTimeCallback([500.0], i -> begin
    i.p[idx] = 50.0
    set_proposed_dt!(i, 0.01)
end)
cb2 = PresetTimeCallback([1000.0], i -> begin
    i.p[idx] = 100.0
    set_proposed_dt!(i, 0.01)
end)
cb3 = PresetTimeCallback([1500.0], i -> begin
    i.p[idx] = 150.0
    set_proposed_dt!(i, 0.01)
end)
cb4 = PresetTimeCallback([2000.0], i -> begin
    i.p[idx] = 0.0
    set_proposed_dt!(i, 0.01)
end)

prob = ODEProblem(osys, [], (0.0, 2500.0))
sol = solve(prob, callback=CallbackSet(cb1, cb2, cb3, cb4))

@unpack M, Y, L = osys
fig = plot(sol, idxs=[Y], xlabel="Time (min)", title="Fig 7.7 (A)", label="β-galactosidase monomer")

lac = function (t)
    if 500 < t < 1000
        50
    elseif 1000 < t < 1500
        100
    elseif 1500 < t < 2000
        150
    else
        0
    end
end

plot!(fig, lac, 0, 2500, label="External lactose (μM)")

fig |> PNG

# ## Fig 7.07 (B)
# Compare the original model and the modified model

rn_mod = @reaction_network begin
    a1 / (1 + RToverK1 * (K2 / (K2 + L))^4), 0 --> M
    δM, M --> 0
    (c1 * M, δY), 0 <--> Y
    mm(Le, kL * 4 * Enz, KML), 0 --> L
    mm(L, 2 * kg * Enz, KMg), L ⇒ 0
    δL, L --> 0
end

#---
setdefaults!(rn_mod, [
    :δM => 0.48,
    :δY => 0.03,
    :δL => 0.02,
    :a1 => 0.29,
    :K2 => 2.92 * 1e6,
    :RToverK1 => 213.2,
    :c1 => 18.8,
    :kL => 6 * 1e4,
    :KML => 680,
    :kg => 3.6 * 1e3,
    :KMg => 7 * 1e5,
    :Le => 0.0,
    :M => 0.01,
    :Y => 0.1,
    :L => 0.0,
    :Enz => 40.0
])

osys_mod = convert(ODESystem, rn_mod; remove_conserved=true)
equations(osys_mod)

#---
prob = SteadyStateProblem(rn, [])
prob_mod = SteadyStateProblem(rn_mod, [])

@unpack Le = prob.f.sys
idx = findfirst(isequal(Le), parameters(prob.f.sys))

@unpack Le = prob_mod.f.sys
idx_mod = findfirst(isequal(Le), parameters(prob_mod.f.sys))

@unpack Y = prob.f.sys
lerange = range(0, 100, 101)

eprob = EnsembleProblem(prob;
    prob_func=(prob, i, repeat) -> begin
        prob.p[idx] = lerange[i]
        prob
    end,
    output_func=(sol, i) -> (sol[Y] / 4, false)
)

eprob_mod = EnsembleProblem(prob_mod;
    prob_func=(prob, i, repeat) -> begin
        prob.p[idx_mod] = lerange[i]
        prob
    end,
    output_func=(sol, i) -> (sol[Y] / 4, false)
)

sim = solve(eprob, DynamicSS(Rodas5()); trajectories=length(lerange))
sim_mod = solve(eprob_mod, DynamicSS(Rodas5()); trajectories=length(lerange))

fig = plot(lerange, sim.u, label="Original")
fig = plot!(fig, lerange, sim_mod.u, label="Modified")
fig = plot!(fig,
    xlabel="External lactose concentration (μM)",
    ylabel="β-galactosidase",
    title="Fig 7.7 (B)"
)

fig |> PNG
