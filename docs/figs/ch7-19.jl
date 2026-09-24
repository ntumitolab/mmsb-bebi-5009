# # Fig 7.19
# Circadian rhythm model
using DifferentialEquations
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

@time "Build system" rn719 = @reaction_network begin
    hillr(PN, vs, ki, n), 0 --> M
    mm(M, vm, km1), M => 0
    ks * M, 0 --> P0
    mm(P0, v1, k1), P0 => P1
    mm(P1, v2, k2), P1 => P0
    mm(P1, v3, k3), P1 => P2
    mm(P2, v4, k4), P2 => P1
    k1, P2 --> PN
    k2, PN --> P2
    mm(P2, vd, kd), P2 => 0
end

#---
alg = FBDF()
up719 = Dict(
    :vs => 0.76, :vm => 0.65, :vd => 0.95, :ks => 0.38, # :kt1 => 1.9, # :kt2 => 1.3,
    :v1 => 3.2, :v2 => 1.58, :v3 => 5.0, :v4 => 2.5,
    :k1 => 1.0, :k2 => 1.0, :k3 => 2.0, :k4 => 2.0, :ki => 1.0, :km1 => 0.5, :kd => 0.2, :n => 4,
    :M => 1.0, :P0 => 1.0, :P1 => 0.0, :P2 => 0.0, :PN => 0.0
)
tspan = (-50.0, 200.0)
@time "Build problem" prob719 = ODEProblem(rn719, up719, tspan)
@time "Solve problem" sol719 = solve(prob719, FBDF())
@unpack M, P0, P1, P2, PN = rn719
totalP = P0 + P1 + P2 + PN
plot(sol719, idxs=[M, PN, totalP], labels=["M" "Nuclear PER" "Total PER"], xlabel="Time", ylabel="Concentration", title="Fig 7.19 (A)")
