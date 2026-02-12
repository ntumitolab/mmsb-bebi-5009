#===
# Fig 3.03

Michaelis-Menten kinetics
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box)

# Enzyme kinetics full model
function model303(; tend = 1.0)
    @parameters k1=30 km1=1 k2=10 ET=1
    @variables S(t)=5.0 ES(t)=0.0 P(t)=0.0 E(t)

    v1 = k1 * S * E - km1 * ES
    v2 = k2 * ES

    eqs = [
        E ~ ET - ES
        D(S) ~ -v1
        D(ES) ~ v1 - v2
        D(P) ~ v2
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end
@time prob = model303()

#---
@time sol = solve(prob, Tsit5())

#---
@unpack S, ES, P, E = prob.f.sys
plot(sol, idxs=[S, ES, E, P], xlabel="Time", ylabel="Concentration", title="Fig. 3.03")

# ## QSSA of ES complex
function model303mm(s0=5; tend=1.0)
    @parameters k1=30 km1=1 k2=10 ET=1 S0=s0
    @variables P(t)=0.0 S(t)
    eqs = [
        S ~ S0 - P
        D(P) ~ (k2 * k1 * S * ET) / (km1 + k2 + k1 * S)
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

@time prob303mm = model303mm()

#---
@time sol303mm = solve(prob303mm, Tsit5())

#---
plot(sol, idxs=[S, P], xlabel="Time", ylabel="Concentration", title="Fig. 3.03 (QSSA)", label=["S (full)" "P (full)"], linestyle=:dash)
plot!(sol303mm, idxs=[S, P], label=["S (QSSA)" "P (QSSA)"], legend=:right)
