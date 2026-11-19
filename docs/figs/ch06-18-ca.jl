# # Fig 6.18
# Model of calcium-induced calcium release in hepatocytes, without Catalyst:
using ComponentArrays
using DiffEqCallbacks
using OrdinaryDiffEq
using Plots
using Startup: hil
Plots.gr(linewidth=1.5)

#---
## We write each flux once as a function of the state and parameters, and
## `cicr_ode!` just balances every species: production minus consumption (the
## same `du = S * v` split Catalyst uses internally).
function cicr_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1p = p.k1 * p.I * u.R                          ## R --> RI
    v1m = p.km1 * u.RI                              ## RI --> R
    v2p = p.k2 * u.C * u.RI                         ## RI --> RIC
    v2m = p.km2 * u.RIC                             ## RIC --> RI
    v3p = p.k3 * u.C * u.RIC                        ## RIC --> RICC
    v3m = p.km3 * u.RICC                            ## RICC --> RIC
    v4 = p.vr * (p.γ0 + p.γ1 * u.RIC) * (p.Cer - u.C)  ## 0 --> C
    v5 = p.p1 * hil(u.C, p.p2, 4)                     ## C => 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.R = -v1p + v1m
    du.RI = v1p - v1m - v2p + v2m
    du.RIC = v2p - v2m - v3p + v3m
    du.RICC = v3p - v3m
    du.C = v4 - v5
    return nothing
end

#---
p618 = (k1=12.0, km1=8.0, k2=15.0, km2=1.65, k3=1.8, km3=0.21, vr=0.185, γ0=0.1, γ1=20.5, p1=8.5, p2=0.065, Cer=8.37)
u618 = ComponentArray(C=0.0, R=1.0, RI=0.0, RIC=0.0, RICC=0.0)
alg = FBDF()

# ## Fig 6.18 (A)
# Constant IP3 level: the discrete parameter `I` is just a component of `p`.
p618a = ComponentArray(; p618..., I=2.0)
tspan = (0.0, 25.0)
@time "Build problem" prob618a = ODEProblem(cicr_ode!, u618, tspan, p618a)
@time "Solve problem" sol618a = solve(prob618a, alg)
plot(sol618a, idxs=[:C, :RIC, :RICC], labels=["C" "RIC" "RICC"], title="Fig 6.18 (A)", xlabel="Time", ylabel="Abundance", legend=:topright)

# ## Fig 6.18 (B)
# Stepped IP3 input, replacing the three `SymbolicDiscreteCallback`s that set
# `I ~ 0.7`, `I ~ 1.2`, and `I ~ 4.0` at `t = 20, 60, 90`. `PresetTimeCallback`
# samples the input signal `I_at` at each jump time and makes the solver end
# its step there, so `I` stays piecewise constant between events.
I_at(t) = t < 20.0 ? 0.0 : t < 60.0 ? 0.7 : t < 90.0 ? 1.2 : 4.0
update_I!(integrator) = integrator.p.I = I_at(integrator.t)
cbs = PresetTimeCallback([20.0, 60.0, 90.0], update_I!)

p618b = ComponentArray(; p618..., I=0.0)
tend = 120.0
@time "Build problem" prob618b = ODEProblem(cicr_ode!, u618, (0.0, tend), p618b)
@time "Solve problem" sol618b = solve(prob618b, alg; callback=cbs)
plot(sol618b, idxs=[:C], label="Calcium", title="Fig 6.18 (B)", xlabel="Time", ylabel="Concentration", legend=:topright)
