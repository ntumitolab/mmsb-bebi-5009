# # Fig 6.3
# Two component pathway
using OrdinaryDiffEq
using Catalyst
using Plots

@time "Build system" rn603 = @reaction_network begin
    @discrete_events begin
        (t == 1.0) => [L => 3.0]
        (t == 3.0) => [L => 0.0]
    end
    (k1 * L, km1), R <--> RL
    k2, RL + P --> RL + Ps
    k3, Ps --> P
end

#---
alg = FBDF()

up603 = Dict(:k1 => 5.0, :km1 => 1.0, :k2 => 6.0, :k3 => 2.0, :L => 0.0, :R => 3.0, :RL => 0.0, :P => 8.0, :Ps => 0.0)
tend = 10.0
@time "Build problem" prob603 = ODEProblem(rn603, up603, (0.0, tend); remove_conserved=true)
@time "Solve problem" sol603 = solve(prob603, alg; tstops=[1.0, 3.0])
plot(sol603, idxs=[:RL, :Ps, :L], labels=["RL" "P*" "L"], title="Fig. 6.3 (A)", xlabel="Time", ylabel="Concentration")
