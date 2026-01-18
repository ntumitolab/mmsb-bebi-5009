#===
# Stochastic simulations

## Gillespie Algorithm
===#
using CairoMakie
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
fig, ax, sp = series(soldirect.t, soldirect.u', labels=["A", "B"]; axis=(;  xlabel="time", ylabel="# of molecules", title="SSA (direct method)"))
axislegend(ax)
fig

# Plot the solution by the first reaction method
fig, ax, sp = series(solfirst.t, solfirst.u', labels=["A", "B"]; axis=(;  xlabel="time", ylabel="# of molecules", title="SSA (first reaction method)"))
axislegend(ax)
fig

# Running 50 simulations
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
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "# of molecules",
    title = "SSA (first method) ensemble"
)
for sol in sols
    lines!(ax, sol.t, sol.u[:, 1], color = (:blue, 0.05))
    lines!(ax, sol.t, sol.u[:, 2], color = (:red, 0.05))
end
lines!(ax, ts, a_avg.(ts), color = :black, linewidth = 3, linestyle = :solid, label = "Average [A]")
lines!(ax, ts, b_avg.(ts), color = :black, linewidth = 3, linestyle = :dash, label = "Average [B]")
axislegend(ax, position = :rt)
fig
