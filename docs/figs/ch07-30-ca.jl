# # Fig. 7.30
# model of synthetic band detector system
using ComponentArrays
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
using Startup: hil, hilr, with_param
Plots.gr(linewidth=1.5)

#---
## The reactions
## ```
## hillr(L, aG, KL, 2),            0 --> G
## bG,                             G --> 0
## hillr(C, aL1, KC, 2),           0 --> L
## mm(R, aL2, KR),                 0 --> L
## bL,                             L --> 0
## mm(R, aC, KR),                  0 --> C
## bC,                             C --> 0
## k1 * (RT - 2R)^2 * A^2,         0 --> R
## k2,                             R --> 0
## ```
function band_detector_ode!(du, u, p, t)
    ## Flux through each reaction:
    w1 = p.aG * hilr(u.L, p.KL, 2)             ## 0 --> G, repressed by L
    w2 = p.bG * u.G                            ## G --> 0
    w3 = p.aL1 * hilr(u.C, p.KC, 2)            ## 0 --> L, repressed by C
    w4 = p.aL2 * hil(u.R, p.KR)                ## 0 --> L, activated by R
    w5 = p.bL * u.L                            ## L --> 0
    w6 = p.aC * hil(u.R, p.KR)                 ## 0 --> C, activated by R
    w7 = p.bC * u.C                            ## C --> 0
    w8 = p.k1 * (p.RT - 2u.R)^2 * p.A^2        ## 0 --> R, driven by A
    w9 = p.k2 * u.R                            ## R --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.G = w1 - w2
    du.L = w3 + w4 - w5
    du.C = w6 - w7
    du.R = w8 - w9
    return nothing
end

#---
p = ComponentArray(
    aG=2.0,      ## muM/min
    bG=0.07,     ## /min
    KL=0.8,      ## muM
    aL1=1.0,     ## muM/min
    aL2=1.0,     ## muM/min
    KC=0.008,    ## muM
    bL=0.02,     ## /min
    aC=1.0,      ## muM/min
    bC=0.07,     ## /min
    k1=0.5,      ## /muM^3 /min
    k2=0.02,     ## /min
    KR=0.01,     ## muM
    RT=0.5,      ## muM
    A=0.1,       ## muM
)
u0 = ComponentArray(G=0.0, L=0.0, C=0.0, R=0.0)

#---
(sols, as) = let N = 1:100
    as = [exp10(-4 + 4 * (i - 1) / length(N)) for i in N]
    @time "Build problem" prob = SteadyStateProblem(band_detector_ode!, u0, p)
    ## Sweeping `A` means rebuilding the whole parameter array: pair-based
    ## `remake(prob; p=[:A => a])` silently no-ops for a `SteadyStateProblem`.
    prob_func = (prob, ctx) -> remake(prob; p = with_param(prob.p, :A, as[ctx.sim_id]))
    eprob = EnsembleProblem(prob; prob_func)
    @time "Solve problem" sols = solve(eprob, DynamicSS(FBDF()); trajectories=length(N))
    (sols, as)
end;

## Each trajectory's `s.u` is a `ComponentArray`, so `s.u[:L]` selects the
## species by name (the solution object itself is not labelled: `s[:L]` would
## error).
luxI = map(s -> s.u[:L], sols.u)
CI = map(s -> s.u[:C], sols.u)
GFP = map(s -> s.u[:G], sols.u)
plot(as, [luxI, CI, GFP], xscale=:log10, xlabel="AHL concentration (μM)", ylabel="Concentration (μM)", title="Fig. 7.30", label=["LuxI" "cI" "GFP"])
