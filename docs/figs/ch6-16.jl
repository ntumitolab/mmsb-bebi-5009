# # Fig 6.16
# Model of apoptosis signalling pathway
using Catalyst
using OrdinaryDiffEq
using Plots

@time "Build system" rn616 = @reaction_network begin
    @discrete_events begin
        (t == 100.0) => [I => 200.0]
        (t == 1200.0) => [I => 0.0]
    end
    (k1, k2), 0 <--> C8
    k3 * (C3s + I), C8 --> C8s
    k4, C8s --> 0
    (k5, k6), C8s + BAR <--> C8sBAR
    (k7, k8), 0 <--> C3
    k9 * C8s, C3 --> C3s
    k10, C3s --> 0
    (k11, k12), C3s + IAP <--> C3sIAP
    (k13, k14), 0 <--> BAR
    (k15, k16 + k17 * C3s), 0 <--> IAP
    k18, C8sBAR --> 0
    k19, C3sIAP --> 0
end

#---
alg = FBDF()
up616 = Dict(
    :C8 => 1.3E5,
    :C8s => 0.0,
    :C3 => 0.21E5,
    :C3s => 0.0,
    :BAR => 0.4E5,
    :IAP => 0.4E5,
    :C8sBAR => 0.0,
    :C3sIAP => 0.0,
    :k1 => 507.0,
    :k2 => 3.9e-3,
    :k3 => 1e-5,
    :k4 => 5.8e-3,
    :k5 => 5e-4,
    :k6 => 0.21,
    :k7 => 81.9,
    :k8 => 3.9e-3,
    :k9 => 5.8e-6,
    :k10 => 5.8e-3,
    :k11 => 5e-4,
    :k12 => 0.21,
    :k13 => 40.0,
    :k14 => 1e-3,
    :k15 => 464.0,
    :k16 => 1.16e-2,
    :k17 => 3e-4,
    :k18 => 1.16e-2,
    :k19 => 1.73e-2,
    :I => 0.0,
)

tend = 1800.0
@time "Build problem" prob616 = ODEProblem(rn616, up616, (0.0, tend); remove_conserved=true)
@time "Solve problem" sol616 = solve(prob616, alg; tstops=[100.0, 1200.0])

@unpack C3s, C8s, I = rn616
plot(sol616, idxs=[C3s, C8s, I*100], title="Fig 6.16", xlabel="Time", ylabel="Concentration", legend=:right)
