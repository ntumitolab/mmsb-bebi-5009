#===
# Stochastic simulations

## Gillespie Algorithm
===#

using StatsBase         ## Weights() and sample()
using Plots
using Interpolations
using Statistics        ## mean()
using Random            ## randexp()
import DisplayAs.PNG    ## Save memory and disk space
Random.seed!(2022)


#===
Stochastic chemical reaction: Gillespie Algorithm (direct method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
===#

function ssa_direct(model, u0::AbstractArray, tend, p, stoich; tstart=zero(tend))
    t = tstart   ## Current time
    ts = [t]     ## Time points
    u = copy(u0) ## Current state
    us = copy(u) ## States over time
    while t < tend
        a = model(u, p, t)               ## propensities
        dt = randexp() / sum(a)          ## Time step for the direct method
        du = sample(stoich, Weights(a))  ## Choose the stoichiometry for the next reaction
        u .+= du   ## Update time
        t += dt    ## Update time
        us = [us u]  ## Append state
        push!(ts, t) ## Append time point
    end
    us = collect(us') ## Transpose to make columns as state variables, rows as observations
    return (t = ts, u = us)
end

#===
Stochastic chemical reaction: Gillespie Algorithm (first reaction method)
Adapted from: Chemical and Biomedical Enginnering Calculations Using Python Ch.4-3
===#

function ssa_first(model, u0, tend, p, stoich; tstart=zero(tend))
    t = tstart   ## Current time
    ts = [t]     ## Time points
    u = copy(u0) ## Current state
    us = copy(u) ## States over time
    while t < tend
        a = model(u, p, t)              ## propensities
        dts = randexp(length(a)) ./ a   ## time scales of all reactions
        i = argmin(dts) ## Choose which reaction to occur
        dt = dts[i]
        du = stoich[i]
        u .+= du  ## Update state
        t += dt   ## Update time
        us = [us u]  ## Append state variable to record
        push!(ts, t) ## Append time point to record
    end
    us = collect(us') ## Make column as variables, rows as observations
    return (t = ts, u = us)
end

#===
Propensity model for this example reaction.
Reaction of A <-> B with rate constants k1 & k2
===#

model(u, p, t) = [p.k1 * u[1],  p.k2 * u[2]]

#---
parameters = (k1=1.0, k2=0.5)

# Stoichiometry for each reaction
stoich=[[-1, 1], [1, -1]]

# Initial conditions (Usually discrete values)
u0 = [200, 0]

# Simulation time
tend = 10.0

# Solve the problem using both direct and first reaction method
soldirect = ssa_direct(model, u0, tend, parameters, stoich)
solfirst = ssa_first(model, u0, tend, parameters, stoich)

# Plot the solution from the direct method
plot(soldirect.t, soldirect.u,
    xlabel="time", ylabel="# of molecules",
    title = "SSA (direct method)", label=["A" "B"])

# Plot the solution by the first reaction method
plot(solfirst.t, solfirst.u,
    xlabel="time", ylabel="# of molecules",
    title = "SSA (1st reaction method)", label=["A" "B"])

# Running 50 simulations
numRuns = 50

sols = map(1:numRuns) do _
    ssa_direct(model, u0, tend, parameters, stoich)
end;

# A interpolation function for each solution
itpsA = map(sols) do sol
    linear_interpolation(sol.t, sol.u[:, 1], extrapolation_bc = Line())
end;

itpsB = map(sols) do sol
    linear_interpolation(sol.t, sol.u[:, 2], extrapolation_bc = Line())
end;

# Calculate average A and B levels in the ensemble.
a_avg(t) = mean(i->i(t), itpsA)
b_avg(t) = mean(i->i(t), itpsB)

# Plot the soluton
fig1 = plot(xlabel="Time", ylabel="# of molecules", title = "SSA (direct method) ensemble")

for sol in sols
    plot!(fig1, sol.t, sol.u, linecolor=[:blue :red], linealpha=0.05, label=false)
end

fig1 |> PNG

# Plot averages
plot!(fig1, a_avg, 0.0, tend, linecolor=:black, linewidth=3, linestyle = :solid, label="Avarage [A]")
plot!(fig1, b_avg, 0.0, tend, linecolor=:black, linewidth=3, linestyle = :dash, label="Avarage [B]")

fig1 |> PNG

#===
## Using Catalyst

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) is a domain-specific language (DSL) package to solve law of mass action problems.
===#

using Catalyst

rn = @reaction_network begin
    k1, A --> B
    k2, B --> A
end

# The system with *integer* state variables belongs to `DiscreteProblem`. A `DiscreteProblem` could be further dispatched into other types of problems, such as `ODEProblem`, `SDEProblem`, and `JumpProblem`.

using DifferentialEquations

params = [:k1 => 1.0, :k2 => 0.5]
u0 = [:A => 200, :B => 0]
tend = 10.0

dprob = DiscreteProblem(rn, u0, (0.0, tend), params)

# In this case, we would like to solve a `JumpProblem` using [Gillespie's Direct stochastic simulation algorithm (SSA)](https://doi.org/10.1016/0021-9991(76)90041-3).

jumpProb = JumpProblem(rn, dprob, Direct())
sol = solve(jumpProb, SSAStepper())

using Plots
plot(sol)

# Parallel ensemble simulation
ensprob = EnsembleProblem(jumpProb)
sim = solve(ensprob, SSAStepper(), EnsembleThreads(); trajectories=50)
plot(sim, alpha=0.1, color=[:blue :red]) |> PNG

summ = EnsembleSummary(sim, 0:0.1:10)
plot(summ,fillalpha=0.5) |> PNG

#===
**See also** the [JumpProcesses.jl docs](https://docs.sciml.ai/JumpProcesses/stable/api/#JumpProcesses.ConstantRateJump) about discrete stochastic examples.
- High-level solutions using `Catalyst.jl` and low-level solutions defining the jumps directly.
- Coupling stochastic discrete jumping and ODEs.
- `RegularJumps` using a more efficient tau-leaping method.
- [More solvers](https://docs.sciml.ai/JumpProcesses/stable/jump_types/) for discrete stochastic simulations.
===#
