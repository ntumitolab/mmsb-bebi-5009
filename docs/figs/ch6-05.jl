# # Fig 6.5
# Model of G-protein signalling pathway
using OrdinaryDiffEq
using SteadyStateDiffEq
using Catalyst
using Plots
Plots.gr(linewidth=1.5)
# ## Fig 6.5 (A)
@time "Build system" rn605 = @reaction_network begin
    @discrete_events begin
        (t == 200.0) => [L => 1e-9]
        (t == 800.0) => [L => 0.0]
    end
    kRL * L, R --> RL
    kRLm, RL --> R
    kGa, RL + G --> RL + Ga + Gbg
    kGd0, Ga --> Gd
    kG1, Gd + Gbg --> G
end

#---
alg = FBDF()
tend = 1200.0
up605 = Dict(:R => 4e3, :RL => 0.0, :G => 1e4, :Ga => 0.0, :Gd => 0.0, :Gbg => 0.0, :kRL => 2e6, :kRLm => 0.01, :kGa => 1e-5, :kGd0 => 0.11, :kG1 => 1.0, :L => 0.0)
@time "Build problem" prob605 = ODEProblem(rn605, up605, (0.0, tend))
@time "Solve problem" sol = solve(prob605, alg; tstops=[200.0, 800.0])
plot(sol, idxs=[:RL, :Ga], labels=["RL" "Ga"], title="Fig. 6.05 (A)", xlabel="Time", ylabel="Concentration")

# ## Fig 6.5 (B)
lrange = range(0, 20 * 1e-9, length=101)
@time "Build system" rn605b = @reaction_network begin
    kRL * L, R --> RL
    kRLm, RL --> R
    kGa, RL + G --> RL + Ga + Gbg
    kGd0, Ga --> Gd
    kG1, Gd + Gbg --> G
end

@time "Build problem" prob605b = SteadyStateProblem(rn605b, up605)

@time "Ensemble simulations" sols605 = map(lrange) do lval
    _p = remake(prob605b; p=[:L => lval])
    sol = solve(_p, DynamicSS(alg))
end

ga = [sol[:Ga] for sol in sols605]
rl = [sol[:RL] for sol in sols605]
plot(lrange .* 1e9, [ga , rl], label=["Ga" "RL"], xlabel="L (nM)", ylabel="Steady-state abundance", title = "Fig. 6.5 (B)")
