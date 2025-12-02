#===
# Stochastic simulations

## Gillespie Algorithm
===#
using Plots
using DataInterpolations: LinearInterpolation
using StatsBase: Weights, sample
using Statistics: mean
using Random
using DisplayAs: PNG    ## For faster rendering
Random.seed!(2024)

#===
## Do-it-yourself

Stochastic chemical reaction: Gillespie Algorithm (direct and first reaction method)
Adapted from: Chemical and Biomedical Engineering Calculations Using Python Ch.4-3
===#
function ssa_alg(model, u0::AbstractVector, tend, p, stoich; tstart=zero(tend), method=:direct)
    t = tstart   ## Current time
    ts = [t]     ## Time points
    u = copy(u0) ## Current state
    us = copy(u') ## States over time
    while t < tend
        a = model(u, p, t)               ## propensities
        if method == :direct
            dt = randexp() / sum(a)          ## Time step for the direct method
            du = sample(stoich, Weights(a))  ## Choose the stoichiometry for the next reaction
        elseif method == :first
            dts = randexp(length(a)) ./ a   ## time scales of all reactions
            i = argmin(dts)                 ## Choose the most recent reaction to occur
            dt = dts[i]
            du = stoich[i]
        else
            error("Method should be either :direct or :first")
        end
        u .+= du   ## Update time
        t += dt    ## Update time
        us = [us; u']  ## Append state
        push!(ts, t)   ## Append time point
    end
    return (t=ts, u=us)
end

#===
Propensity model for this example reaction.
Reaction of A <-> B with rate constants k1 & k2
===#

model(u, p, t) = [p.k1 * u[1], p.k2 * u[2]]

#---
parameters = (k1=1.0, k2=0.5)

# Stoichiometry for each reaction
stoich = [[-1, 1], [1, -1]]

# Initial conditions (Usually discrete values)
u0 = [200, 0]

# Simulation time
tend = 10.0

# Solve the problem using both direct and first reaction method
@time soldirect = ssa_alg(model, u0, tend, parameters, stoich; method=:direct)
@time solfirst = ssa_alg(model, u0, tend, parameters, stoich; method=:first)

# Plot the solution from the direct method
plot(soldirect.t, soldirect.u,
    xlabel="time", ylabel="# of molecules",
    title="SSA (direct method)", label=["A" "B"]) |> PNG

# Plot the solution by the first reaction method
plot(solfirst.t, solfirst.u,
    xlabel="time", ylabel="# of molecules",
    title="SSA (1st reaction method)", label=["A" "B"]) |> PNG

# Running 50 simulations
numRuns = 50

@time sols = map(1:numRuns) do _
    ssa_alg(model, u0, tend, parameters, stoich; method=:direct)
end;

# Average values and interpolation
ts = range(0, tend, 101)
a_avg(t) = mean(sols) do sol
    A = LinearInterpolation(sol.u[:, 1], sol.t; cache_parameters=true)
    A(t)
end

b_avg(t) = mean(sols) do sol
    A = LinearInterpolation(sol.u[:, 2], sol.t; cache_parameters=true)
    A(t)
end

# Plot the solution
fig1 = plot(xlabel="Time", ylabel="# of molecules", title="SSA (first method) ensemble")

for sol in sols
    plot!(fig1, sol.t, sol.u, linecolor=[:blue :red], linealpha=0.05, label=false)
end

fig1 |> PNG

# Plot averages
plot!(fig1, a_avg, 0.0, tend, linecolor=:black, linewidth=3, linestyle=:solid, label="Average [A]") |> PNG
plot!(fig1, b_avg, 0.0, tend, linecolor=:black, linewidth=3, linestyle=:dash, label="Average [B]") |> PNG
fig1 |> PNG

#===
## Using Catalyst (recommended)

https://github.com/SciML/Catalyst.jl is a domain-specific language (DSL) package to simulate chemical reaction networks.
===#
using Catalyst
using JumpProcesses
using Plots

two_state_model = @reaction_network begin
    k1, A --> B
    k2, B --> A
end

# State variables are *integers*

params = [:k1 => 1.0, :k2 => 0.5]
u0 = [:A => 200, :B => 0]
tspan = (0.0, 10.0)
prob = DiscreteProblem(two_state_model, u0, tspan, params)
jump_prob = JumpProblem(two_state_model, prob, Direct())

# In this case, we would like to solve a `JumpProblem` using [Gillespie's Direct stochastic simulation algorithm (SSA)](https://doi.org/10.1016/0021-9991(76)90041-3).
@time sol = solve(jump_prob, SSAStepper())
plot(sol) |> PNG

# Parallel ensemble simulation
ensprob = EnsembleProblem(jump_prob)
@time sim = solve(ensprob, SSAStepper(), EnsembleThreads(); trajectories=50)

#---
plot(sim, alpha=0.1, color=[:blue :red]) |> PNG

#---
summ = EnsembleSummary(sim, 0:0.1:10)
plot(summ, fillalpha=0.5) |> PNG

# ### SIR model
sir_model = @reaction_network begin
    beta, S + I --> 2I
    gamma, I --> R
end

p = (:beta => 0.1 / 1000, :gamma => 0.01)
u0 = [:S => 990, :I => 10, :R => 0]
tspan = (0.0, 250.0)
prob = DiscreteProblem(sir_model, u0, tspan, p)
jump_prob = JumpProblem(sir_model, prob, Direct())
#---
@time sol = solve(jump_prob, SSAStepper())

#---
plot(sol) |> PNG

#===
**See also** the [JumpProcesses.jl docs](https://docs.sciml.ai/JumpProcesses/stable/) about discrete stochastic algorithm examples.
===#

# ## Grid simulation
# https://docs.sciml.ai/Catalyst/stable/spatial_modelling/lattice_reaction_systems/
# https://docs.sciml.ai/JumpProcesses/stable/tutorials/spatial/
using Catalyst
using JumpProcesses
using Plots

sir_model = @reaction_network begin
    beta, S + I --> 2I
    gamma, I --> R
end

dS = @transport_reaction D S
dI = @transport_reaction D I
lattice = CartesianGrid((3,3))
lrs = LatticeReactionSystem(sir_model, [dS, dI], lattice)

s0 = ones(Int, 3, 3) .* 110
i0 = zeros(Int, 3, 3)
i0[1, 1] = 10
r0 = zeros(Int, 3, 3)
u0 = [:S => s0, :I => i0, :R => r0]
ps = [:beta => 0.1 / 100, :gamma => 0.03, :D => 1.0]
tspan = (0.0, 250.0)
prob = DiscreteProblem(lrs, u0, tspan, ps)
jump_prob = JumpProblem(lrs, prob, NSM())

@time sol = solve(jump_prob, SSAStepper())
#---
lat_getu(sol, :S, lrs)
#---
lat_getu(sol, :I, lrs)
#---
lat_getu(sol, :R, lrs)

#---
t = sol.t
s = sum.(lat_getu(sol, :S, lrs))
i = sum.(lat_getu(sol, :I, lrs))
r = sum.(lat_getu(sol, :R, lrs))
plot(t, [s i r], label=["S" "I" "R"], lw=1.5) |> PNG
