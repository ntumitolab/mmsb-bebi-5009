# methionine model

using ComponentArrays
using Plots
using DifferentialEquations
using UnPack

ps = ComponentArray(
    v_MATI_max=561., K_MATI_m=41., K_MATI_i=50.,
    v_MATIII_max=22870., K_MATIII_m2=21.1,
    v_MET_max=4544., A_over_K_MET_m2 = 0.1,
    v_GNMT_max=10600., K_GNMT_m=4500., K_GNMT_i=20.,
    alpha_d=1333., K_AHC=0.1, Adenosine=1., Met=48.5
)

function metmodel!(D, u, p, t)
    @unpack AdoMet, AdoHcy = u
    K_MATIII_m1 = 20000 / (1 + 5.7 * (AdoMet / (AdoMet + 600)^2))
    K_MET_m1 = 10 * (1 + AdoHcy / 4)
    @unpack K_AHC, Adenosine = p
    Hcy = AdoHcy * K_AHC / Adenosine


end
