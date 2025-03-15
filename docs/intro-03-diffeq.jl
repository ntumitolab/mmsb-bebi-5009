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
using Plots
Plots.default(linewidth=2)

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
plot(sol)

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
using Plots
Plots.default(linewidth=2)

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
plot(sol, label=["S" "I" "R"], legend=:right)

#===
## Using ModelingToolkit.jl (recommended)

[ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/) is a high-level package for symbolic-numeric modeling and simulation in the Julia ecosystem.
===#
using ModelingToolkit
using OrdinaryDiffEq
using Plots
Plots.default(linewidth=2)

# ### Exponential decay model

@independent_variables t    ## Time
@parameters λ               ## Decaying rate constant
@variables C(t)             ## Time and concentration
D = Differential(t)         ## Differential operator

# Define an ODE with equations
eqs = [D(C) ~ -λ * C]
@mtkbuild decaySys = ODESystem(eqs, t)

#---
u0 = [C => 1.0]
p = [λ => 1.0]
tspan = (0.0, 2.0)

prob = ODEProblem(decaySys, u0, tspan, p)
sol = solve(prob)
plot(sol)

# ### SIR model
using OrdinaryDiffEq
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

@independent_variables t
@parameters β γ
@variables s(t) i(t) r(t)
D = Differential(t)

eqs = [
    D(s) ~ -β * s * i,
    D(i) ~ β * s * i - γ * i,
    D(r) ~ γ * i
]

@mtkbuild sirSys = ODESystem(eqs, t)

#---
p = [β => 1.0, γ => 0.3]
u0 = [s => 0.99, i => 0.01, r => 0.00]
tspan = (0.0, 20.0)

prob = ODEProblem(sirSys, u0, tspan, p)
sol = solve(prob)

plot(sol)

#===
### Lorenz system

The Lorenz system is a system of ordinary differential equations having chaotic solutions for certain parameter values and initial conditions. ([Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system))

$$
\begin{align}
  \frac{dx}{dt} &= \sigma(y-x) \\
  \frac{dy}{dt} &= x(\rho - z) -y \\
  \frac{dz}{dt} &= xy - \beta z
\end{align}
$$
===#
using OrdinaryDiffEq
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

# Building the model is wrapped in a function
function build_lorentz(; name)
    @parameters begin
        σ = 10.0
        ρ = 28.0
        β = 8 / 3
    end

    @independent_variables t
    @variables begin
        x(t) = 1.0    ## Independent variable (time)
        y(t) = 0.0    ## Independent variable (time)
        z(t) = 0.0    ## Independent variable (time)
    end

    D = Differential(t)

    eqs = [
        D(x) ~ σ * (y - x),
        D(y) ~ x * (ρ - z) - y,
        D(z) ~ x * y - β * z
    ]

    sys = ODESystem(eqs, t; name)
    return sys
end

#---
tspan = (0.0, 100.0)
@mtkbuild sys = build_lorentz()
prob = ODEProblem(sys, [], tspan, [])
sol = solve(prob)

# x-y-z time-series
plot(sol)

# y-t plot
plot(sol, idxs=[sys.y])

# `idxs=(sys.x, sys.y, sys.z)` makes a phase plot.
plot(sol, idxs=(sys.x, sys.y, sys.z), label=false, size=(600, 600))

# ## Saving simulation results
using DataFrames
using CSV

df = DataFrame(sol)
CSV.write("lorenz.csv", df)
rm("lorenz.csv")

#===
## Catalyst.jl

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) is a domain-specific language (DSL) package to simulate chemical reaction networks.
===#
using OrdinaryDiffEq
using Catalyst
using Plots
Plots.default(linewidth=2)

# ### Exponential decay model
decay_rn = @reaction_network begin
    λ, C --> 0
end

#---
p = [:λ => 1.0]
u0 = [:C => 1.0]
tspan = (0.0, 2.0)

prob = ODEProblem(decay_rn, u0, tspan, p)
sol = solve(prob)

#---
plot(sol, title="Exponential Decay")

# ### SIR model
sir_rn = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
end

# Extract the symbols for later use
@unpack β, γ, S, I, R = sir_rn

p = [β => 1.0, γ => 0.3]
u0 = [S => 0.99, I => 0.01, R => 0.00]
tspan = (0.0, 20.0)

prob = ODEProblem(sir_rn, u0, tspan, p)
sol = solve(prob)

#---
plot(sol, legend=:right, title="SIR Model")
