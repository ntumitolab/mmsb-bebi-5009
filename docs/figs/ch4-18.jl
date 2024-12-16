#===
# Fig 4.18

Continuation diagram

See also [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl)
===#
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---
function model418!(D, u, p, t)
    a, b = u
    k1, k2, k3, k4, k5, n = p
    D[1] = k1 / (1 + b^n) - (k3 + k5) * a
    D[2] = k2 + k5 * a - k4 * b
end

function ainf(k1)
    ps = (k1 = k1, k2 = 5., k3 = 5., k4 = 5., k5 = 2., n = 4.)
    u0 = [0., 0.]
    prob = SteadyStateProblem(model418!, u0, ps)
    sol = solve(prob, DynamicSS(Rodas5()))
    return sol.u[1]
end

#---
plot(
    ainf, 0., 1000.,
    title = "Fig 4.18",
    xlabel = "K1" , ylabel= "Steady state [A]",
    legend=nothing, ylim=(0, 4), xlim=(0, 1000)
)
