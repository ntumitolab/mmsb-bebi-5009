#===
# Fig 6.18

Model of calcium-induced calcium release in hepatocytes
===#
using Catalyst
using ModelingToolkit
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---
rn = @reaction_network begin
    (k1 * I, km1), R <--> RI
    (k2 * C, km2), RI <--> RIC
    (k3 * C, km3), RIC <--> RICC
    vr * (γ0 + γ1 * RIC) * (Cer - C), 0 --> C
    hill(C, p1, p2, 4), C ⇒ 0
end

#---
setdefaults!(rn, [
    :γ0 => 0.1,
    :γ1 => 20.5,
    :p1 => 8.5,
    :p2 => 0.065,
    :k1 => 12,
    :k2 => 15,
    :k3 => 1.8,
    :km1 => 8,
    :km2 => 1.65,
    :km3 => 0.21,
    :Cer => 8.37,
    :vr => 0.185,
    :I => 0,
    :C => 0,
    :R => 1,
    :RI => 0,
    :RIC => 0,
    :RICC => 0
])

osys = convert(ODESystem, rn; remove_conserved = true) |> structural_simplify
equations(osys)

# ## Fig 6.18 (A)
@unpack I = osys
prob = ODEProblem(osys, [], (0., 25.), [I => 2.0])
sol = solve(prob)

@unpack C, RIC, RICC = osys
plot(sol, idxs=[C, RIC, RICC], title="Fig 6.18 (A)", xlabel="Time", ylabel="Abundance", legend=:topright)

# ## Fig 6.18 (B)
cb1 = PresetTimeCallback([20.0], i -> begin i.ps[I] = 0.7; set_proposed_dt!(i, 0.01) end)
cb2 = PresetTimeCallback([60.0], i -> begin i.ps[I] = 1.2; set_proposed_dt!(i, 0.01) end)
cb3 = PresetTimeCallback([90.0], i -> begin i.ps[I] = 4.0; set_proposed_dt!(i, 0.01) end)
cbs = CallbackSet(cb1, cb2, cb3)

prob = ODEProblem(osys, [], (0., 120.))
sol = solve(prob, callback=cbs)

plot(sol, idxs=[osys.C], title="Fig 6.18 (B)", xlabel="Time", ylabel="Ca concentration", legend=false, ylim=(0, 2.5))
