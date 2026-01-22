# # Figure 7.30
# model of synthetic band detector system
using OrdinaryDiffEq
using SteadyStateDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

#---
function model730!(D, u, p, t)
    hil(x, k) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @unpack aG, bG, KL, aL1, aL2, KC, bL, aC, bC, k1, k2, KR, RT, A = p
    @unpack G, L, C, R = u
    D.G = -bG * G + aG * hil(KL, L, 2)
    D.L = -bL * L + aL1 * hil(KC, C, 2) + aL2 * hil(R, KR)
    D.C = -bC * C + aC * hil(R, KR)
    D.R = -k2 * R + k1 * (RT - 2 * R)^2 * A^2
    nothing
end

p0_730 = ComponentArray(
    aG=2.0, # muM/min
    bG=0.07, # /min
    KL=0.8, # muM
    aL1=1, # muM/min
    aL2=1, # muM/min
    KC=0.008, # muM
    bL=0.02, # /min
    aC=1, # muM/min
    bC=0.07, # /min
    k1=0.5, # /muM^3 /min
    k2=0.02, # /min
    KR=0.01, # muM
    RT=0.5, # muM
    A=0.1,
)

u0 = ComponentArray(
    G=0.0,
    L=0.0,
    C=0.0,
    R=0.0,
)

#---
N = 1:100
as = [exp10(-4 + 4*(i-1)/length(N)) for i in N]
prob = SteadyStateProblem(model730!, u0, p0_730)
prob_func = (prob, i, repeat) -> remake(prob, p=ComponentArray(p0_730; A=as[i]))
eprob = EnsembleProblem(prob; prob_func)

@time sol = solve(eprob, DynamicSS(KenCarp47()); trajectories=length(N))

#---
luxI = map(s -> s.u.L, sol)
CI = map(s -> s.u.C, sol)
GFP = map(s -> s.u.G, sol)

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "External AHL concentration (μM)",
    ylabel = "Concentration (μM)",
    title = "Fig 7.30",
    xscale = log10,
)
lines!(ax, as, GFP, label = "GFP")
lines!(ax, as, CI, label = "cI")
lines!(ax, as, luxI, label = "LuxI")
axislegend(ax, position = :lc)
fig
