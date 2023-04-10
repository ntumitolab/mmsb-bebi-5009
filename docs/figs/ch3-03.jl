#===
# Fig 3.03

Michaelis-Menten kinetics
===#

using DifferentialEquations
using Catalyst
using Plots
Plots.default(linewidth=2)

#---
rn303 = @reaction_network begin
    (k1, km1), S + E <--> ES
    k2, ES --> E + P
end

#---
u0 = [:S=>5., :ES=>0., :E=>1., :P=>0.]
ps = [:k1 => 30., :km1 => 1., :k2 => 10.]
tend = 1.0

#---
prob = ODEProblem(rn303, u0, tend, ps)
sol = solve(prob)

#---
plot(sol, xlabel="Time (AU)", ylabel="Concentration (AU)", legend=:right, title="Fig 3.03")

#---
rn303mm = @reaction_network begin
    mm(S, k2 * ET, (km1 + k2) / k1), S ⇒ P ## using \Rightarrow
end

#---
u0 = [:S=>5., :P=>0.]
ps = [:k1 => 30., :km1 => 1., :k2 => 10., :ET=>1.]
tend = 1.0
probmm = ODEProblem(rn303mm, u0, tend, ps)
solmm = solve(probmm)

#---
@unpack S, P = rn303
fig = plot(sol, idxs=[S, P], line=(:dash), label=["S (full)" "P (full)"])
plot!(fig, solmm, idxs=[S, P], label=["S (MM)" "P (MM)"])
plot!(fig, title="Fig. 3.03",
    xlabel="Time (AU)", ylabel="Concentration (AU)",
    xlims=(0., tend), ylims=(0., 5.), legend=:right)

# ## Runtime information
import Pkg
Pkg.status()

#---
import InteractiveUtils
InteractiveUtils.versioninfo()
