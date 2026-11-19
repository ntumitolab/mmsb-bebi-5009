# # Fig 6.16
# Model of apoptosis signalling pathway (caspase cascade with positive feedback
# and drug inhibition)
using ComponentArrays
using DiffEqCallbacks
using OrdinaryDiffEq
using Plots
using Startup: mass_action
Plots.gr(linewidth=1.5)

#---
## The reactions, exactly as written in the Catalyst version:
## ```
## (k1, k2),                0 <--> C8
## k3 * (C3s + I),          C8 --> C8s
## k4,                      C8s --> 0
## (k5, k6),                C8s + BAR <--> C8sBAR
## (k7, k8),                0 <--> C3
## k9 * C8s,                C3 --> C3s
## k10,                     C3s --> 0
## (k11, k12),              C3s + IAP <--> C3sIAP
## (k13, k14),              0 <--> BAR
## (k15, k16 + k17 * C3s),  0 <--> IAP
## k18,                     C8sBAR --> 0
## k19,                     C3sIAP --> 0
## ```
function apoptosis_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1p = mass_action(p.k1)                          ## 0 --> C8
    v1m = mass_action(p.k2, u.C8)                    ## C8 --> 0
    v2 = mass_action(p.k3 * (u.C3s + p.I), u.C8)     ## C8 --> C8s, activated by C3s and drug I
    v3 = mass_action(p.k4, u.C8s)                    ## C8s --> 0
    v4p = mass_action(p.k5, u.C8s, u.BAR)            ## C8s + BAR --> C8sBAR
    v4m = mass_action(p.k6, u.C8sBAR)                ## C8sBAR --> C8s + BAR
    v5p = mass_action(p.k7)                          ## 0 --> C3
    v5m = mass_action(p.k8, u.C3)                    ## C3 --> 0
    v6 = mass_action(p.k9 * u.C8s, u.C3)             ## C3 --> C3s, activated by C8s
    v7 = mass_action(p.k10, u.C3s)                   ## C3s --> 0
    v8p = mass_action(p.k11, u.C3s, u.IAP)           ## C3s + IAP --> C3sIAP
    v8m = mass_action(p.k12, u.C3sIAP)               ## C3sIAP --> C3s + IAP
    v9p = mass_action(p.k13)                         ## 0 --> BAR
    v9m = mass_action(p.k14, u.BAR)                  ## BAR --> 0
    v10p = mass_action(p.k15)                        ## 0 --> IAP
    v10m = mass_action(p.k16 + p.k17 * u.C3s, u.IAP) ## IAP --> 0, degraded faster in the presence of C3s
    v11 = mass_action(p.k18, u.C8sBAR)               ## C8sBAR --> 0
    v12 = mass_action(p.k19, u.C3sIAP)               ## C3sIAP --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.C8 = v1p - v1m - v2
    du.C8s = v2 - v3 - v4p + v4m
    du.BAR = v9p - v9m - v4p + v4m
    du.C8sBAR = v4p - v4m - v11
    du.C3 = v5p - v5m - v6
    du.C3s = v6 - v7 - v8p + v8m
    du.IAP = v10p - v10m - v8p + v8m
    du.C3sIAP = v8p - v8m - v12
    return nothing
end

I_signal(t) = 200.0 * (100.0 <= t < 1200.0)
affect_I!(integrator) = (integrator.p.I = I_signal(integrator.t))
cbs = PresetTimeCallback([100.0, 1200.0], affect_I!)

#---
p616 = ComponentArray(
    k1=507.0, k2=3.9e-3, k3=1e-5, k4=5.8e-3, k5=5e-4, k6=0.21,
    k7=81.9, k8=3.9e-3, k9=5.8e-6, k10=5.8e-3, k11=5e-4, k12=0.21,
    k13=40.0, k14=1e-3, k15=464.0, k16=1.16e-2, k17=3e-4, k18=1.16e-2, k19=1.73e-2,
    I=0.0,
)
u616 = ComponentArray(
    C8=1.3E5, C8s=0.0, BAR=0.4E5, C8sBAR=0.0,
    C3=0.21E5, C3s=0.0, IAP=0.4E5, C3sIAP=0.0,
)

tend = 1800.0
@time "Build problem" prob616 = ODEProblem(apoptosis_ode!, u616, (0.0, tend), p616)
@time "Solve problem" sol616 = solve(prob616, FBDF(), callback=cbs)

#---
plot(sol616, idxs=[:C3s, :C8s], title="Fig 6.16", xlabel="Time", ylabel="Concentration", legend=:right)
plot!(sol616.t, 100 .* I_signal.(sol616.t), st=:steppost, label="I * 100")
