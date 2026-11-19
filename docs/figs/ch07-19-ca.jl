# # Fig 7.19 A
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: hil, hilr
Plots.gr(linewidth=1.5)

#---
## The reactions
## ```
## hillr(PN, vs, ki, n), 0 --> M
## mm(M, vm, km1),       M => 0
## ks * M,               0 --> P0
## mm(P0, v1, k1),       P0 => P1
## mm(P1, v2, k2),       P1 => P0
## mm(P1, v3, k3),       P1 => P2
## mm(P2, v4, k4),       P2 => P1
## k1,                   P2 --> PN
## k2,                   PN --> P2
## mm(P2, vd, kd),       P2 => 0
## ```
function circadian_ode!(du, u, p, t)
    ## Flux through each reaction:
    w1 = p.vs * hilr(u.PN, p.ki, p.n)  ## 0 --> M
    w2 = p.vm * hil(u.M, p.km1)        ## M => 0
    w3 = p.ks * u.M                    ## 0 --> P0
    w4 = p.v1 * hil(u.P0, p.k1)        ## P0 => P1
    w5 = p.v2 * hil(u.P1, p.k2)        ## P1 => P0
    w6 = p.v3 * hil(u.P1, p.k3)        ## P1 => P2
    w7 = p.v4 * hil(u.P2, p.k4)        ## P2 => P1
    w8 = p.k1 * u.P2                   ## P2 --> PN
    w9 = p.k2 * u.PN                   ## PN --> P2
    w10 = p.vd * hil(u.P2, p.kd)       ## P2 => 0

    ## Species balances: each flux enters with its stoichiometric coefficient.
    ## `w8` leaves `P2` and enters `PN` (nuclear import), and `w9` comes back
    ## the other way, so the two species trade protein without any loss.
    du.M = w1 - w2
    du.P0 = w3 + w5 - w4
    du.P1 = w4 + w7 - w5 - w6
    du.P2 = w6 + w9 - w7 - w8 - w10
    du.PN = w8 - w9
    return nothing
end

#---
## The rate constants of the Catalyst model, in its `up719` dictionary: the
## parameters of the fluxes live in the parameter array `p`, the five species
## in the state array `u0`. `n = 4` is the Hill coefficient of the feedback.
p = ComponentArray(
    vs=0.76, vm=0.65, vd=0.95, ks=0.38,
    v1=3.2, v2=1.58, v3=5.0, v4=2.5,
    k1=1.0, k2=1.0, k3=2.0, k4=2.0, ki=1.0, km1=0.5, kd=0.2, n=4.0,
)
u0 = ComponentArray(M=1.0, P0=1.0, P1=0.0, P2=0.0, PN=0.0)

alg = FBDF()
tspan = (-50.0, 200.0)
@time "Build problem" prob719 = ODEProblem(circadian_ode!, u0, tspan, p)
@time "Solve problem" sol719 = solve(prob719, alg)

#---
ts = range(tspan[1], tspan[2]; length=500)
totalP = sol719(ts; idxs=:P0) + sol719(ts; idxs=:P1) + sol719(ts; idxs=:P2) + sol719(ts; idxs=:PN)

## Component names work as indices into the solution too: `idxs=[:M, :PN]`
## selects the species by name, without knowing their order in `u`.
plot(sol719, idxs=[:M, :PN], labels=["M" "Nuclear PER"], xlabel="Time", ylabel="Concentration", title="Fig 7.19 (A)")
plot!(ts, totalP, labels=["Total PER"])
