# methionine model

using LabelledArrays
using Plots
using DifferentialEquations
using UnPack

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

ps = (
    v_MATI_max=561., K_MATI_m=41., K_MATI_i=50.,
    v_MATIII_max=22870., K_MATIII_m2=21.1,
    v_MET_max=4544., A_over_K_MET_m2 = 0.1,
    v_GNMT_max=10600., K_GNMT_m=4500., K_GNMT_i=20.,
    alpha_d=1333., K_AHC=0.1, Adenosine=1., Met=48.5
)

function metmodel(u, p, t)
    @unpack AdoMet, AdoHcy = u
    K_MATIII_m1 = 20000 / (1 + 5.7 * (AdoMet / (AdoMet + 600)^2))
    K_MET_m1 = 10 * (1 + AdoHcy / 4)
    @unpack K_AHC, Adenosine = p
    Hcy = AdoHcy * K_AHC / Adenosine
    @unpack v_MATI_max, Met, K_MATI_m, K_MATI_i = p
    v_MATI = v_MATI_max * hil(Met, K_MATI_m) * hil(K_MATI_i, AdoMet)
    @unpack v_MATIII_max, Met, K_MATIII_m2 = p
    v_MATIII = v_MATIII_max * hil(Met, K_MATIII_m1 * hil(K_MATIII_m2, Met))
    @unpack v_GNMT_max, K_GNMT_m, K_GNMT_i = p
    v_GNMT = v_GNMT_max * hil(AdoMet, K_GNMT_m, 2.3) * hil(K_GNMT_i, AdoHcy)
end

function metmodel!(D, u, p, t)
    D.AdoMet, D.AdoHcy = metmodel(u, p, t)
    return nothing
end
