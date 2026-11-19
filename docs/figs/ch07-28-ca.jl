# # Fig 7.28
# model of synthetic pulse generating system
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: hil
Plots.gr(linewidth=1.5)

#---
## Hill-type activation of GFP production by LuxR:AHL (`R`) and cI (`C`),
## copied from the Catalyst version.
function vg728(R, C, KR, KC, aG)
    fR = R / KR
    fC = C / KC
    return aG * (fR / (1 + fR + fC^2 + fR * fC^2))
end

## Michaelis-Menten kinetics come from `Startup`: Catalyst's `mm(x, v, K)` is
## `v * hil(x, K)`.

#---
## The reactions
## ```
## vg728(R, C, KR, KC, aG),        0 --> G
## bG,                             G --> 0
## mm(R, aC, KR),                  0 --> C
## bC,                             C --> 0
## k1 * (RT - 2R)^2 * A^2,         0 --> R
## k2,                             R --> 0
## ```
function pulse_ode!(du, u, p, t)
    ## Flux through each reaction:
    w1 = vg728(u.R, u.C, p.KR, p.KC, p.aG)        ## 0 --> G
    w2 = p.bG * u.G                               ## G --> 0
    w3 = p.aC * hil(u.R, p.KR)                    ## 0 --> C
    w4 = p.bC * u.C                               ## C --> 0
    w5 = p.k1 * (p.RT - 2 * u.R)^2 * p.A^2        ## 0 --> R
    w6 = p.k2 * u.R                               ## R --> 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.G = w1 - w2
    du.C = w3 - w4
    du.R = w5 - w6
    return nothing
end

#---
## Same parameter values as the Catalyst version.
p = ComponentArray(
    aG=80.0,   ## muM/min
    bG=0.07,   ## /min
    KC=0.008,  ## muM
    aC=0.5,    ## muM/min
    bC=0.3,    ## /min
    k1=0.5,    ## /muM^3 /min
    k2=0.02,   ## /min
    KR=0.02,   ## muM
    RT=0.5,    ## muM
    A=10.0,    ## muM
)
u0 = ComponentArray(G=0.0, C=0.0, R=0.0)

alg = FBDF()
tend = 50.0
@time "Build problem" prob728 = ODEProblem(pulse_ode!, u0, (0.0, tend), p)
@time "Solve problem" sol728 = solve(prob728, alg)

#---
## Component names work as indices into the solution too: `idxs=[:G, :C, :R]`
## selects the species by name, without knowing their order in `u`.
plot(sol728, idxs=[:G, :C, :R], xlabel="Time (min)", ylabel="Concentration (μM)", title="Fig 7.28", label=["GFP" "cI" "LuxR:AHL complex"])
