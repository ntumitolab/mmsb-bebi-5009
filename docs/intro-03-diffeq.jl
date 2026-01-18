#===
# Solving differential equations in Julia

**Standard procedures**

- Define a model function representing the right-hand-side (RHS) of the system.
  - Out-of-place form: `f(u, p, t)` where `u` is the state variable(s), `p` is the parameter(s), and `t` is the independent variable (usually time). The output is the right hand side (RHS) of the differential equation system.
  - In-place form: `f!(du, u, p, t)`, where the output is saved to `du`. The rest is the same as the out of place form. The in-place form has potential performance benefits since it allocates less than the out-of-place (`f(u, p, t)`) counterpart.
  - Using ModelingToolkit.jl : define equations and build an ODE system.
- Initial conditions (`u0`) for the state variable(s).
- (Optional) define parameter(s) `p`.
- Define a problem (e.g. `ODEProblem`) using the modeling function (`f`), initial conditions (`u0`), simulation time span (`tspan == (tstart, tend)`), and parameter(s) `p`.
- Solve the problem by calling `solve(prob)`.

## Solve ODEs using OrdinaryDiffEq.jl

Documentation: <https://docs.sciml.ai/DiffEqDocs/stable/>

### Single variable: Exponential decay model

The concentration of a decaying nuclear isotope could be described as an exponential decay:

$$
\frac{d}{dt}C(t) = - \lambda C(t)
$$

**State variable**
- $C(t)$: The concentration of a decaying nuclear isotope.

**Parameter**
- $\lambda$: The rate constant of decay. The half-life $t_{\frac{1}{2}} = \frac{ln2}{\lambda}$

===#
using OrdinaryDiffEq
using CairoMakie

# The model function is the 3-argument out-of-place form, `f(u, p, t)`.
decay(u, p, t) = p * u

p = -1.0            ## Rate of exponential decay
u0 = 1.0            ## Initial condition
tspan = (0.0, 2.0)  ## Start time and end time

prob = ODEProblem(decay, u0, tspan, p)
sol = solve(prob)

# Solution at time t=1.0 (with interpolation)
sol(1.0)

# Time points
sol.t

# Solutions at corresponding time points
sol.u

# Visualize the solution
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Exponential Decay"
)
lines!(ax, sol, label = "C(t)")
axislegend(ax, position = :rt)
fig

#===
### Three variables: The SIR model

The SIR model describes the spreading of an contagious disease can be described by the [SIR model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology):

$$
\begin{align}
\frac{d}{dt}S(t) &= - \beta S(t)I(t)  \\
\frac{d}{dt}I(t) &= \beta S(t)I(t)  - \gamma I(t)  \\
\frac{d}{dt}R(t) &= \gamma I(t)
\end{align}
$$

**State variables**

- $S(t)$ : the fraction of susceptible people
- $I(t)$ : the fraction of infectious people
- $R(t)$ : the fraction of recovered (or removed) people

**Parameters**

- $\beta$ : the rate of infection when susceptible and infectious people meet
- $\gamma$ : the rate of recovery of infectious people

===#
using OrdinaryDiffEq
using CairoMakie

# SIR model (in-place form can save array allocations and thus faster)
function sir!(du, u, p, t)
    s, i, r = u
    β, γ = p
    v1 = β * s * i
    v2 = γ * i
    du[1] = -v1
    du[2] = v1 - v2
    du[3] = v2
    return nothing
end

#---
p = (β=1.0, γ=0.3)
u0 = [0.99, 0.01, 0.00]
tspan = (0.0, 20.0)
prob = ODEProblem(sir!, u0, tspan, p)
sol = solve(prob)

# Visualize the solution
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, 0..20, t-> sol(t)[1], label="S")
lines!(ax, 0..20, t-> sol(t)[2], label="I")
lines!(ax, 0..20, t-> sol(t)[3], label="R")
axislegend(ax, position = :rc)
fig

# ## Saving simulation results
using DataFrames
using CSV

df = DataFrame(sol)
CSV.write("sir.csv", df)
rm("sir.csv")
