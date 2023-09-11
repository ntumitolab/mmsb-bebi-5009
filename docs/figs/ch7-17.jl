#===
# Fig 7.17

Goodwin oscillator model: https://en.wikipedia.org/wiki/Goodwin_model_(biology)
===#

using Catalyst
using ModelingToolkit
using DifferentialEquations
using Plots
Plots.default(linewidth=2)
import DisplayAs.SVG

#---
rn = @reaction_network begin
    (a / (k^n + Z^n), b), 0 <--> X
    (α * X, β), 0 <--> Y
    (γ * Y, δ), 0 <--> Z
end

#---
ps = [
    :a => 360,
    :k => 1.368,
    :b => 1,
    :α => 1,
    :β => 0.6,
    :γ => 1,
    :δ => 0.8,
    :n => 12
]

u0 = zeros(3)
tend = 35.0

#---
prob = ODEProblem(rn, u0, tend, ps)

#---
sol = solve(prob)

#---
fig = plot(sol, title="Fig 7.17 (A)", xlabel="Time", ylabel="Concentration")
fig |> SVG
#---
fig = plot(sol, idxs=(1, 2, 3), title="Fig 7.17 (B)", legend=false, size=(600, 600))

fig |> SVG
