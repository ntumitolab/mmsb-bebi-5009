# # Chapter 6
# ## Fig 6.3
# Two component pathway
using OrdinaryDiffEq
using SteadyStateDiffEq
using Catalyst
using Plots

# Model
@time "Build system" rn603 = @reaction_network begin
    @variables L(t)
    @equations begin
        L ~ 3 * (1 < t) * (t < 3)
    end
    (k1, km1), L + R <--> RL
    k2, RL + P --> RL + Ps
    k3, Ps --> P
end

u0map603 = Dict(:R => 3.0, :RL => 0.0, :P => 8.0, :Ps => 0.0)
psmap603 = Dict(:k1 => 5.0, :km1 => 1.0, :k2 => 6.0, :k3 => 2.0)
tend = 10.0
osys = ode_model(rn603; remove_conserved=true)
@time "Build problem" prob603 = ODEProblem(rn603, u0map603, (0.0, tend), psmap603; remove_conserved=true, mtkcompile=true)
@time
