#===
# Fig 7.7

model of lac operon in E. coli
===#
using ComponentArrays: ComponentArray
using SimpleUnPack
using DiffEqCallbacks
using OrdinaryDiffEq
using SteadyStateDiffEq
using CairoMakie

#---
function model707!(D, u, p, t)
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @unpack Le, a1, RToverK1, K2, δM, δY, δL, c1, kL, KML, kg, KMg = p
    @unpack M, Y, L = u
    D.M = a1 / (1 + RToverK1 * (K2 / (K2 + L))^4) - δM * M
    D.Y = c1 * M - δY * Y
    D.L = kL * Y * hil(Le, KML) - 2kg * (Y/4) * hil(L, KMg) - δL * L
    nothing
end

#---
ps707 = ComponentArray(
    δM = 0.48,
    δY = 0.03,
    δL = 0.02,
    a1 = 0.29,
    K2 = 2.92e6,
    RToverK1 = 213.2,
    c1 = 18.8,
    kL = 6e4,
    KML = 680.0,
    kg = 3.6e3,
    KMg = 7e5,
    Le = 0.0,
)

ics707 = ComponentArray(
    M = 0.01,
    Y = 0.1,
    L = 0.0,
)

# Events
affect_le1!(integrator) = integrator.p.Le = 50
affect_le2!(integrator) = integrator.p.Le = 100
affect_le3!(integrator) = integrator.p.Le = 150
affect_le4!(integrator) = integrator.p.Le = 0
event_le1 = PresetTimeCallback([500.0], affect_le1!)
event_le2 = PresetTimeCallback([1000.0], affect_le2!)
event_le3 = PresetTimeCallback([1500.0], affect_le3!)
event_le4 = PresetTimeCallback([2000.0], affect_le4!)
cbs = CallbackSet(event_le1, event_le2, event_le3, event_le4)

# ## Fig 7.07 (A)
tend = 2500.0
prob707a = ODEProblem(model707!, ics707, tend, ps707)
@time sol = solve(prob707a, TRBDF2(), callback=cbs)

#---
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (min)", ylabel="Concentration", title="Fig 7.7 (A)\nLac operon model in E. coli")
lines!(ax, 0..tend, t -> sol(t).Y, label="β-galactosidase monomer")
lines!(ax, 0..tend, t -> 50 * (500<=t<1000) + 100 * (1000<=t<1500) + 150 * (1500<=t<2000), label="External lactose")
axislegend(ax, position=:lt)
fig

# ## Fig 7.07 (B)
# Compare the original model and the modified model
function model707b!(D, u, p, t)
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @unpack Le, a1, RToverK1, K2, δM, δY, δL, c1, kL, KML, kg, KMg, Enz = p
    @unpack M, Y, L = u
    D.M = a1 / (1 + RToverK1 * (K2 / (K2 + L))^4) - δM * M
    D.Y = c1 * M - δY * Y
    D.L = 4kL * Enz * hil(Le, KML) - 2kg * Enz * hil(L, KMg) - δL * L
    nothing
end

#---
ps707b = ComponentArray(ps707; Enz=40.0)
prob707a = SteadyStateProblem(model707!, ics707, ps707)
prob707b = SteadyStateProblem(model707b!, ics707, ps707b)

lerange = range(0, 100, 101)

eprob = EnsembleProblem(prob707a;
    prob_func=(prob, i, repeat) -> remake(prob, p=ComponentArray(ps707; Le=lerange[i])),
    output_func=(sol, i) -> (sol.u.Y / 4, false)
)

eprob_mod = EnsembleProblem(prob707b;
    prob_func=(prob, i, repeat) -> remake(prob, p=ComponentArray(ps707b; Le=lerange[i])),
    output_func=(sol, i) -> (sol.u.Y / 4, false)
)

alg = DynamicSS(TRBDF2())
@time sim = solve(eprob, alg; trajectories=length(lerange))
@time sim_mod = solve(eprob_mod, alg; trajectories=length(lerange))

fig = Figure()
ax = Axis(fig[1, 1], xlabel="External lactose concentration (μM)", ylabel="β-galactosidase", title="Fig 7.7 (B)\nLac operon model in E. coli")
lines!(ax, lerange, sim.u, label="Original model")
lines!(ax, lerange, sim_mod.u, label="Modified model")
axislegend(ax, position=:lt)
fig
