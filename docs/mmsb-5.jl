# # Chapter 5
# Methionine model

import DisplayAs.PNG
using DifferentialEquations
using LabelledArrays
using UnPack
using Plots
Plots.default(linewidth=2)

# Convenience functions

hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

# ## Methionine model

function metmodel(u, p, t)
    AdoMet, AdoHcy = u

    @unpack K_AHC, Adenosine = p
    Hcy = AdoHcy * K_AHC / Adenosine
    @unpack v_MATI_max, Met, K_MATI_m, K_MATI_i = p
    v_MATI = v_MATI_max * hil(Met * hil(K_MATI_i, AdoMet), K_MATI_m)
    @unpack v_MATIII_max, Met, K_MATIII_m2 = p
    K_MATIII_m1 = 20000 / (1 + 5.7 * hil(AdoMet, 600)^2)
    v_MATIII = v_MATIII_max * hil(Met, K_MATIII_m1 * hil(K_MATIII_m2, Met))
    @unpack v_GNMT_max, K_GNMT_m, K_GNMT_i = p
    v_GNMT = v_GNMT_max * hil(AdoMet, K_GNMT_m, 2.3) * hil(K_GNMT_i, AdoHcy)
    @unpack v_MET_max, A_over_K_MET_m2 = p
    K_MET_m1 = 10 + 2.5 * AdoHcy
    v_MET = v_MET_max * hil(AdoMet, K_MET_m1) * hil(A_over_K_MET_m2)
    @unpack alpha_d = p
    v_D = alpha_d * Hcy

    dAdoMetdt = (v_MATI + v_MATIII) - (v_GNMT + v_MET)
    dAdoHcydt = (v_GNMT + v_MET - v_D) * hil(Adenosine, K_AHC)
    return (dAdoMetdt, dAdoHcydt)
end

function metmodel!(D, u, p, t)
    D[1], D[2] = metmodel(u, p, t)
    return nothing
end

#---

ps = (
    v_MATI_max=561., K_MATI_m=41., K_MATI_i=50.,
    v_MATIII_max=22870., K_MATIII_m2=21.1,
    v_MET_max=4544., A_over_K_MET_m2 = 0.1,
    v_GNMT_max=10600., K_GNMT_m=4500., K_GNMT_i=20.,
    alpha_d=1333., K_AHC=0.1, Adenosine=1., Met=48.5
)

tend = 5.
u0 = LVector(AdoMet=10., AdoHcy=10.)

prob = ODEProblem(metmodel!, u0, tend, ps)
sol = solve(prob)

# ## Figure 5.10

fig510 = plot(sol, title="Figure 5.10", xlabel="Time (hr)", ylabel="Concentration (μM)", xlims=(0, 1), legend=:right)

fig510 |> PNG

# ## Figure 5.11A

rx = range(0, 1200, 101)
ry = range(0, 6, 101)

∂A = (x, y) -> metmodel((x, y), ps, 0)[1]
∂B = (x, y) -> metmodel((x, y), ps, 0)[2]

fig511a = plot(title="Figure 5.11A")
contour!(fig511a, rx, ry, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig511a, identity, 0, 0, line=(:black), label="AdoMet nullcline")
contour!(fig511a, rx, ry, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig511a, identity, 0, 0, line=(:black, :dash), label="AdoHcy nullcline")

tend = 15.
u0s = (
    [500.,1.5],
    [900.,2.5],
    [1100.,3.5],
    [400.,5.0],
    [800.,5.5],
    [1000.,5.75],
    [300.,1],
    [700.,2],
    [200.,5],
    [600.,5.25]
)

for u0 in u0s
    local prob = ODEProblem(metmodel!, u0, tend, ps)
    local sol = solve(prob)
    plot!(fig511a, sol, idxs=(1, 2), label=false, alpha=0.5)
end

plot!(fig511a, xlims=(0, 1200), ylims=(0, 6), xlabel="AdoMet (μM)", ylabel="AdoHcy (μM)", legend=:bottomright)

fig511a |> PNG

# ## Figure 5.11B

# Methionine level increased.

ps2 = merge(ps, (;Met=51))

rx = range(0, 1200, 101)
ry = range(0, 6, 101)

∂A = (x, y) -> metmodel((x, y), ps2, 0)[1]
∂B = (x, y) -> metmodel((x, y), ps2, 0)[2]

fig511b = plot(title="Figure 5.11B")
contour!(fig511b, rx, ry, ∂A, levels=[0], cbar=false, line=(:black))
plot!(fig511b, identity, 0, 0, line=(:black), label="AdoMet nullcline")
contour!(fig511b, rx, ry, ∂B, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig511b, identity, 0, 0, line=(:black, :dash), label="AdoHcy nullcline")

tend = 15.
u0s = (
    [420,1.5],
    [820,2.5],
    [1120,3.5],
    [520,5.0],
    [620,2],
    [240,4.5],
    [720,5.5]
)

for u0 in u0s
    plot!(fig511b, ODEProblem(metmodel!, u0, tend, ps) |> solve, idxs=(1, 2), label=false, alpha=0.5)
end

plot!(fig511b, xlims=(0, 1200), ylims=(0, 6), xlabel="AdoMet (μM)", ylabel="AdoHcy (μM)", legend=:bottomright)

fig511b |> PNG

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()
