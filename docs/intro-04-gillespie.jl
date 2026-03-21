#===
# Stochastic simulations

## Gillespie Algorithm
===#
using Plots
using DataInterpolations: LinearInterpolation
using StatsBase: Weights, sample
using Statistics: mean
using Random
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
plot(soldirect.t, soldirect.u, labels=["A" "B"], xlabel="time", ylabel="# of molecules", title="SSA (direct method)")

# Plot the solution by the first reaction method
plot(solfirst.t, solfirst.u, labels=["A" "B"], xlabel="time", ylabel="# of molecules", title="SSA (first reaction method)")

# Running 50 simulations serially
# TODO: parallelize this loop using `Threads.@threads` or `Distributed` for better performance
@time sols = map(1:50) do _
    ssa_alg(model, u0, tend, parameters, stoich; method=:direct)
end;

# Average values and interpolation
ts = range(0, tend, 101)
a_interp = [LinearInterpolation(sol.u[:, 1], sol.t; cache_parameters=true) for sol in sols]
b_interp = [LinearInterpolation(sol.u[:, 2], sol.t; cache_parameters=true) for sol in sols]
a_avg(t) = mean(i-> i(t), a_interp)
b_avg(t) = mean(i-> i(t), b_interp)

# Plot the solution
plot(ts, a_avg.(ts), label="Average [A]", line=(:black, :solid))
plot!(ts, b_avg.(ts), label="Average [B]", line=(:black, :dash))
for sol in sols
    plot!(sol.t, sol.u[:, 1], line = (:blue, 0.05), label=false)
    plot!(sol.t, sol.u[:, 2], line = (:red, 0.05), label=false)
end
plot!(xlabel="time", ylabel="# of molecules", title="SSA (first method) ensemble")
