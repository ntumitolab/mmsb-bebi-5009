#===
# Fig 3.03

Michaelis-Menten kinetics
===#
using OrdinaryDiffEq
using ModelingToolkit
using Catalyst
using Plots
Plots.default(linewidth=2)

# Reaction network
rn = @reaction_network begin
    (k1, km1), S + E <--> ES
    k2, ES --> E + P
end

#---
setdefaults!(rn, [
    :S=>5.,
    :ES=>0.,
    :E=>1.,
    :P=>0.,
    :k1 => 30.,
    :km1 => 1.,
    :k2 => 10.,
])

osys = convert(ODESystem, rn; remove_conserved = true) |> structural_simplify

#---
observed(osys)

#---
tend = 1.0

#---
prob = ODEProblem(osys, [], tend)
sol = solve(prob)

#---
@unpack S, ES, E, P = osys
plot(sol, idxs=[S, ES, E, P], xlabel="Time (AU)", ylabel="Concentration (AU)", legend=:right, title="Fig 3.03")

#---
rn303mm = @reaction_network begin
    mm(S, k2 * ET, (km1 + k2) / k1), S ⇒ P ## using \Rightarrow
end

#---
setdefaults!(rn303mm, [
    :S=>5.,
    :ET=>1.,
    :P=>0.,
    :k1 => 30.,
    :km1 => 1.,
    :k2 => 10.,
])

osysmm = convert(ODESystem, rn303mm; remove_conserved = true) |> structural_simplify

#---
tend = 1.0
probmm = ODEProblem(osysmm, [], tend)
solmm = solve(probmm)

#---
@unpack S, P = osys
fig = plot(sol, idxs=[S, P], line=(:dash), label=["S (full)" "P (full)"])
plot!(fig, solmm, idxs=[S, P], label=["S (MM)" "P (MM)"])
plot!(fig, title="Fig. 3.03",
    xlabel="Time (AU)", ylabel="Concentration (AU)",
    xlims=(0., tend), ylims=(0., 5.), legend=:right
)
