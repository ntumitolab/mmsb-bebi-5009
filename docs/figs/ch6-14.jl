#===
# Fig 6.14

Model of E. coli chemotaxis signalling pathway
===#

using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---

rn = @reaction_network begin
    mm(Am, k1 * BP, KM1), Am ⇒ A
    mm(AmL, k2 * BP, KM2), AmL ⇒ AL
    km1 * R , A ⇒ Am
    km2 * R , AL ⇒ AmL
    (k3 * L, km3), Am <--> AmL
    (k4 * L, km4), A <--> AL
    (k5 * Am, km5), B <--> BP
end

#---

setdefaults!(rn, [
    :Am => 0.0360,
    :AmL => 1.5593,
    :A => 0.0595,
    :AL => 0.3504,
    :B => 0.7356,
    :BP => 0.2644,
    :k1 => 200,
    :k2 => 1,
    :k3 => 1,
    :km1 => 1,
    :km2 => 1,
    :km3 => 1,
    :k5 => 0.05,
    :km5 => 0.005,
    :k4 => 1,
    :km4 => 1,
    :KM1 => 1,
    :KM2 => 1,
    :R => 1,
    :L => 20
])

osys = convert(ODESystem, rn; remove_conserved = true)

#---
observed(osys)

#---
equations(osys)

#---
@unpack L = osys
idx = findfirst(isequal(L), parameters(osys))

#---
cb1 = PresetTimeCallback([10.0], i -> begin i.p[idx] = 40; set_proposed_dt!(i, 0.01) end)
cb2 = PresetTimeCallback([30.0], i -> begin i.p[idx] = 80; set_proposed_dt!(i, 0.01) end)
cbs = CallbackSet(cb1, cb2)

prob = ODEProblem(osys, [], (0., 50.))

#---
sol = solve(prob, callback=cbs)

@unpack Am = osys
fig = plot(sol, idxs=[Am], title="Fig 6.14", xlabel="Time", ylabel="Active CheA ([Am])", ylims=(0.01, 0.04), legend=false)
