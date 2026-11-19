# # Fig 2.11-14
# Intro to model reduction of ODE metabolic networks.
using ComponentArrays
using OrdinaryDiffEq
using Plots
Plots.gr(linewidth=1.5)

#---
## Each reaction has a mass-action flux: its rate constant times the
## product of the concentrations of its reactants.
mass_action(k) = k          ## zero-order flux (constant production)
mass_action(k, x) = k * x   ## first-order flux

# ## Fig 2.11
## The full model
## ```
## k0,  0 --> A
## k1,  A --> B
## km1, B --> A
## k2,  B --> 0
## ```
## and its balances, one flux at a time:
function full_ode!(du, u, p, t)
    v0 = mass_action(p.k0)        ## 0 --> A
    v1 = mass_action(p.k1, u.A)   ## A --> B
    v2 = mass_action(p.km1, u.B)  ## B --> A
    v3 = mass_action(p.k2, u.B)   ## B --> 0

    du.A = v0 - v1 + v2
    du.B = v1 - v2 - v3
    return nothing
end

#---
p211 = ComponentArray(k0=0.0, k1=9.0, km1=12.0, k2=2.0)
u211 = ComponentArray(A=0.0, B=10.0)
tend = 3.0
@time "Build problem" prob211 = ODEProblem(full_ode!, u211, (0.0, tend), p211)
@time "Solve problem" sol211 = solve(prob211, Tsit5())
plot(sol211, idxs=[:A, :B], labels=["A" "B"], title="Fig 2.11", xlabel="Time", ylabel="Concentration")

# ## Fig 2.12
## Rapid equilibrium assumption. The fluxes of `A <--> B` (rate constants 9
## and 12) are much larger than the fluxes feeding and draining the pool, so
## the pair stays close to its fast equilibrium `k1 * A ~ km1 * B` at all
## times. With the conserved total `C = A + B`, the two species track the
## pool through the equilibrium fractions
## ```
## A ~ C * km1 / (km1 + k1)
## B ~ C * k1  / (km1 + k1)
## ```
## and the reduced network acts on `C` alone:
## ```
## k0, 0 --> C
## k2 * k1 / (km1 + k1), C --> 0
## ```
frac_a(p) = p.km1 / (p.km1 + p.k1)  ## fraction of the pool held as A
frac_b(p) = p.k1 / (p.km1 + p.k1)   ## fraction of the pool held as B
a_of_c(c, p) = frac_a(p) * c        ## algebraic reconstruction of A
b_of_c(c, p) = frac_b(p) * c        ## algebraic reconstruction of B

function rapid_eq_ode!(du, u, p, t)
    du.C = mass_action(p.k0) - mass_action(p.k2 * frac_b(p), u.C)
    return nothing
end

#---
p212 = ComponentArray(k0=0.0, k1=9.0, km1=12.0, k2=2.0)
u212 = ComponentArray(C=10.0)  ## the pool A + B, i.e. A(0) + B(0)
@time "Build problem" prob212 = ODEProblem(rapid_eq_ode!, u212, (0.0, tend), p212)
@time "Solve problem" sol212 = solve(prob212, Tsit5())

## `A` and `B` come back from `C` algebraically, at any times we ask for:
ts = range(0.0, tend; length=200)
c212 = sol212(ts, idxs=:C)
plot(sol211, idxs=[:A, :B], linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(ts, [a_of_c.(c212, Ref(p212)) b_of_c.(c212, Ref(p212))],
     linestyle=:solid, label=["A (Reduced model)" "B (Reduced model)"])
plot!(title="Fig 2.12", xlabel="Time", ylabel="Concentration")

# ## Fig 2.13
## When another set of parameters breaks the rapid equilibrium assumption:
## the input flux `k0` is now comparable to the fluxes that equilibrate
## `A <--> B`, and the initial condition `(A, B) = (8, 4)` is far from the
## fast equilibrium of the pool `C = 12`.
p213 = ComponentArray(k0=9.0, k1=21.0, km1=12.0, k2=2.0)

@time "Build problem" prob213full = ODEProblem(full_ode!, ComponentArray(A=8.0, B=4.0), (0.0, tend), p213)
@time "Build problem" prob213re = ODEProblem(rapid_eq_ode!, ComponentArray(C=12.0), (0.0, tend), p213)
@time sol213full = solve(prob213full, Tsit5())
@time sol213re = solve(prob213re, Tsit5())

c213 = sol213re(ts, idxs=:C)
plot(sol213full, idxs=[:A, :B], linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(ts, [a_of_c.(c213, Ref(p213)) b_of_c.(c213, Ref(p213))],
     linestyle=:solid, label=["A (Reduced model)" "B (Reduced model)"])
plot!(title="Fig 2.13", xlabel="Time", ylabel="Concentration")

# ## Fig 2.14
## Quasi-steady state assumption on species `A` (`d(A) ~ 0`). Setting the
## balance of `A` to zero,
## ```
## 0 ~ k0 - (k1 * A - km1 * B)
## ```
## solves for the algebraic concentration `A ~ (k0 + km1 * B) / k1`.
## Substituting it into the balance for `B` gives the reduced model directly:
## ```
## d(B) = (k1 * A - km1 * B) - k2 * B
##      = (k0 + km1 * B) - km1 * B - k2 * B
##      = k0 - k2 * B
## ```
## (Catalyst writes the same substitution as `A ~ (k0 + km1 * B) / k1` in
## `@equations`; ModelingToolkit as the algebraic equation `0 ~ d(A)`.)
a_qss(b, p) = (p.k0 + p.km1 * b) / p.k1  ## algebraic concentration of A

function qss_ode!(du, u, p, t)
    du.B = mass_action(p.k0) - mass_action(p.k2, u.B)
    return nothing
end

#---
p214 = ComponentArray(k0=9.0, k1=21.0, km1=12.0, k2=2.0)
## Start on the slow manifold, at the QSS split of the pool `A + B = 12`
## from Fig 2.13: solving `a_qss(b, p) + b = 12` for `b` gives
b0_qss(c, p) = (p.k1 * c - p.k0) / (p.km1 + p.k1)
@time "Build problem" prob214 = ODEProblem(qss_ode!, ComponentArray(B=b0_qss(12.0, p214)), (0.0, tend), p214)
@time sol214 = solve(prob214, Tsit5())

plot(sol213full, idxs=[:A, :B], linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(ts, [a_qss.(sol214(ts, idxs=:B), Ref(p214)) sol214(ts, idxs=:B)],
     linestyle=:solid, label=["A (QSS model)" "B (QSS model)"])
plot!(title="Fig 2.14", xlabel="Time", ylabel="Concentration")
