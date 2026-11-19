# # Fig 4.18
# Continuation diagram.
# See also [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl) for a more robust continuation diagram.
using ComponentArrays
using OrdinaryDiffEq
using SteadyStateDiffEq
using Plots
Plots.gr(linewidth=1.5)

#---
## The reactions
## ```
## (k1 / (1 + B^n), k3), 0 <--> A
## k5, A --> B
## (k2, k4), 0 <--> B
## ```
function switch_ode!(du, u, p, t)
    v1 = p.k1 / (1 + u.B^p.n)   ## 0 --> A, inhibited by B
    v3 = p.k3 * u.A             ## A --> 0
    v5 = p.k5 * u.A             ## A --> B
    v2 = p.k2                   ## 0 --> B
    v4 = p.k4 * u.B             ## B --> 0

    du.A = v1 - v3 - v5
    du.B = v2 + v5 - v4
    return nothing
end

#---
p = ComponentArray(k1=0.0, k2=5.0, k3=5.0, k4=5.0, k5=2.0, n=4.0)
u0 = ComponentArray(A=0.0, B=0.0)

@time "Build problem" prob = SteadyStateProblem(switch_ode!, u0, p)

#---
with_k1(p, k1) = ComponentArray(merge(NamedTuple(p), (k1 = k1,)))
k1s = 0.0:10.0:1000.0
@time sols = [solve(remake(prob; p=with_k1(p, k1)), DynamicSS(FBDF())).u[:A] for k1 in k1s]
plot(k1s, sols, xlabel="Parameter k1", ylabel="Steady state [A]", title="Fig. 4.18", label=false)
