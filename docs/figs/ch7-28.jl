# # Fig 7.28
# model of synthetic pulse generating system
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

#---
function model728!(D, u, p, t)
    @unpack aG, bG, KC, aC, bC, k1, k2, KR, RT, A = p
    @unpack G, C, R = u
    fR = R / KR
    fC = C / KC
    Gddt = -bG * G + aG * (fR / (1 + fR + fC^2 + fR * fC^2))
    Cddt = -bC * C + aC * R / (KR + R)
    Rddt = -k2 * R + k1 * (RT - 2 * R)^2 * A^2
    D.G = Gddt
    D.C = Cddt
    D.R = Rddt
    nothing
end

#---
ps728 = ComponentArray(
    aG = 80, # muM/min
    bG = 0.07, # /min
    KC= 0.008, # muM
    aC= 0.5, # muM/min
    bC = 0.3, # /min
    k1 = 0.5, # /muM^3 /min
    k2 = 0.02, # /min
    KR= 0.02, # muM
    RT = 0.5, # muM
    A = 10.0, # muM
)

u0728 = ComponentArray(
    G = 0.0,
    C = 0.0,
    R = 0.0,
)

tend = 50.0
prob728 = ODEProblem(model728!, u0728, (0.0, tend), ps728)

#---
@time sol728 = solve(prob728, Tsit5())

#---
plt728 = plot(sol728, labels = ["GFP" "cI" "LuxR:AHL complex"], xlabel="Time (min)", ylabel="Concentration (μM)", title="Fig 7.28", legend=:right)
