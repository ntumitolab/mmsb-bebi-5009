# # Fig 6.14
# Model of E. coli chemotaxis signalling pathway (Barkai–Leibler)
using ComponentArrays
using DiffEqCallbacks
using OrdinaryDiffEq
using Plots
Plots.gr(linewidth=1.5)

#---
## The reactions, exactly as written in the Catalyst version:
## ```
## mm(Am, k1 * BP, KM1),  Am  => A
## mm(AmL, k2 * BP, KM2), AmL => AL
## km1 * R,               A   => Am
## km2 * R,               AL  => AmL
## (k3 * L, km3),         Am  <--> AmL
## (k4 * L, km4),         A   <--> AL
## (k5 * Am, km5),        B   <--> BP
## ```
mm(x, v, k) = v * x / (x + k)  ## Michaelis–Menten factor, as in Catalyst

function chemotaxis_ode!(du, u, p, t)
    ## Flux through each reaction:
    v1 = mm(u.Am, p.k1 * u.BP, p.KM1)  ## Am  => A,   v1 = k1 * BP * Am / (KM1 + Am)
    v2 = mm(u.AmL, p.k2 * u.BP, p.KM2) ## AmL => AL,  v2 = k2 * BP * AmL / (KM2 + AmL)
    v3 = p.km1 * p.R                   ## A   => Am,  zeroth order in A
    v4 = p.km2 * p.R                   ## AL  => AmL, zeroth order in AL
    v5 = p.k3 * p.L * u.Am             ## Am  --> AmL
    v6 = p.km3 * u.AmL                 ## AmL --> Am
    v7 = p.k4 * p.L * u.A              ## A   --> AL
    v8 = p.km4 * u.AL                  ## AL  --> A
    v9 = p.k5 * u.Am * u.B             ## B   --> BP
    v10 = p.km5 * u.BP                 ## BP  --> B

    ## Species balances: each flux enters with its stoichiometric coefficient.
    du.Am = -v1 + v3 - v5 + v6
    du.AmL = -v2 + v4 + v5 - v6
    du.A = v1 - v3 - v7 + v8
    du.AL = v2 - v4 + v7 - v8
    du.B = -v9 + v10
    du.BP = v9 - v10
    return nothing
end

ligand_step(t) = (t < 20.0 ? 40.0 : 80.0)
affect_ligand!(integrator) = (integrator.p.L = ligand_step(integrator.t))
cbs = PresetTimeCallback([10.0, 30.0], affect_ligand!)

#---
p614 = ComponentArray(
    k1=200.0, k2=1.0, k3=1.0, k4=1.0, k5=0.05,
    km1=1.0, km2=1.0, km3=1.0, km4=1.0, km5=0.005,
    KM1=1.0, KM2=1.0, R=1.0, L=20.0,
)
u614 = ComponentArray(Am=0.0360, A=0.0595, AmL=1.5593, AL=0.3504, B=0.7356, BP=0.2644)

tend = 50.0
@time "Build problem" prob614 = ODEProblem(chemotaxis_ode!, u614, (0.0, tend), p614)
@time "Solve problem" sol614 = solve(prob614, FBDF(), callback=cbs)

#---
plot(sol614, idxs=[:Am], title="Fig 6.14", xlabel="Time", ylabel="Concentration")
