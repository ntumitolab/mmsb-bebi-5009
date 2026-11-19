# # Fig 3.03
# Michaelis-Menten kinetics
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: hil
Plots.gr(linewidth=1.5)

#---
# ## Enzyme kinetics full model
# The reaction network
# ```
# k1,  S + E --> ES
# km1, ES    --> S + E
# k2,  ES    --> P + E
# ```
function enzyme_ode!(du, u, p, t)
    ## Flux through each reaction:
    E = p.ET - u.ES          ## free enzyme, from the conservation ET = E + ES
    v1 = p.k1 * u.S * E      ## k1,  S + E --> ES
    v2 = p.km1 * u.ES        ## km1, ES    --> S + E
    v3 = p.k2 * u.ES         ## k2,  ES    --> P + E

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.S = -v1 + v2          ## S is consumed by v1, restored by v2
    du.ES = v1 - v2 - v3     ## ES is made by v1, broken down by v2 and v3
    du.P = v3                ## P is made by v3
    return nothing
end

## Algebraic reconstruction of the eliminated species E.
e_of(es, p) = p.ET - es

#---
p303 = ComponentArray(k1=30.0, km1=1.0, k2=10.0, ET=1.0)
u303 = ComponentArray(S=5.0, ES=0.0, P=0.0)
tend = 1.0
@time "Build problem" prob303 = ODEProblem(enzyme_ode!, u303, (0.0, tend), p303)
@time "Solve problem" sol303 = solve(prob303, Tsit5())

## Component names work as indices into the solution too, and the eliminated
## `E` comes back from the conservation law, so the plot shows all four
## species by name.
ts = range(0.0, tend; length=200)
s_full = vec(sol303(ts; idxs=:S))
es_full = vec(sol303(ts; idxs=:ES))
p_full = vec(sol303(ts; idxs=:P))
e_full = e_of.(es_full, Ref(p303))  ## E = ET - ES
plot(ts, [s_full es_full p_full e_full],
     labels=["S" "ES" "P" "E"], title="Fig. 3.03 (Full model)", xlabel="Time", ylabel="Concentration")

#---
# ## QSSA of the ES complex
# If the complex equilibrates fast compared to substrate turnover, it stays in
# *quasi-steady state*, `d(ES) ~ 0`. Setting its balance to zero and using
# `E = ET - ES`,
# ```
# 0 ~ k1 * S * (ET - ES) - (km1 + k2) * ES
# ```
# solves for `ES ~ ET * S / (Km + S)` with the **Michaelis constant**
# `Km = (km1 + k2) / k1`, and the product-forming flux becomes the
# Michaelis-Menten rate `v = k2 * ES = mm(S, k2 * ET, Km)`. The whole network
# collapses to a single reaction
# ```
# mm(S, k2 * ET, (km1 + k2) / k1), S => P
# ```
# in Catalyst notation, where `=>` means the rate expression is used as-is and
# `mm(X, v, K) = v * X / (X + K)`.
Km(p) = (p.km1 + p.k2) / p.k1            ## Michaelis constant

function qssa_ode!(du, u, p, t)
    ## Flux through the single reaction, used as-is (a `=>` reaction):
    v = (p.k2 * p.ET) * hil(u.S, Km(p))  ## mm(S, k2 * ET, (km1 + k2) / k1), S => P

    du.S = -v                            ## S --> P
    du.P = v                             ## P is made at the same rate
    return nothing
end

#---
p303mm = ComponentArray(k1=30.0, km1=1.0, k2=10.0, ET=1.0)
u303mm = ComponentArray(S=5.0, P=0.0)
@time "Build problem" prob303mm = ODEProblem(qssa_ode!, u303mm, (0.0, tend), p303mm)
@time "Solve problem" sol303mm = solve(prob303mm, Tsit5())

## The reduced model (solid) overlaid on the full model (dashed).
plot(sol303mm, idxs=[:S, :P], labels=["S (QSSA)" "P (QSSA)"], title="Fig. 3.03 (QSSA vs Full)", xlabel="Time", ylabel="Concentration")
plot!(sol303, idxs=[:S, :P], labels=["S (full)" "P (full)"], linestyle=:dash, legend=:right)
