#===

# Solving differential equations in Julia

## Standard procedures

- Define a model function representing the right-hand-side (RHS) of the sysstem.
  - Out-of-place form: `f(u, p, t)` where `u` is the state variable(s), `p` is the parameter(s), and `t` is the independent variable (usually time). The output is the right hand side (RHS) of the differential equation system.
  - In-place form: `f!(du, u, p, t)`, where the output is saved to `du`. The rest is the same as the out of place form. The in-place form has potential performance benefits since it allocates less than the out-of-place (`f(u, p, t)`) counterpart.
- Initial conditions (`u0`) for the state variable(s).
- (Optional) define parameter(s) `p`.
- Define a problem (e.g. `ODEProblem`) using the modeling function (`f`), initial conditions (`u0`), simulation time span (`tspan == (tstart, tend)`), and parameter(s) `p`.
- Solve the problem by calling `solve(prob)`.

## Solve ODEs using DifferentialEquations.jl

Documentation: <https://docs.sciml.ai/DiffEqDocs/stable/>

===#

using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#===

### Exponential decay model

The concentration of a decaying nuclear isotope could be described as an exponential decay:

$$
\frac{d}{dt}C(t) = - \lambda C(t)
$$

**State variable**
- $C(t)$: The concentration of a decaying nuclear isotope.

**Parameter**
- $\lambda$: The rate constant of decay. The half-life $t_{\frac{1}{2}} = \frac{ln2}{\lambda}$

===#

# Model function, in the out-of-place form `f(u, p, t)`

expdecay(u, p, t) = p * u

p = -1.0            ## Rate of exponential decay
u0 = 1.0            ## Initial condition
tspan = (0.0, 2.0)  ## Start time and end time

# Define a problem
prob = ODEProblem(expdecay, u0, tspan, p)

# Solve the problem
sol = solve(prob)

# Visualize the solution
plot(sol, legend=:right)

# Solution at time t=1.0 (with interpolation)
sol(1.0)

# Time points
sol.t

# Solutions at respective time points
sol.u

#===

### SIR model

This 3-variable model describes the spreading of an contagious disease can be described by the [SIR model](https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model):

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

using DifferentialEquations
using Plots
Plots.default(linewidth=2)

# SIR model (in-place form)
function sir!(D, u, p ,t)
	s, i, r = u
	β, γ = p
	v1 = β * s * i
	v2 = γ * i
    D[1] = -v1
    D[2] = v1 - v2
    D[3] = v2
	return nothing
end

#---
p = (β = 1.0, γ = 0.3)
u0 = [0.99, 0.01, 0.00]
tspan = (0.0, 20.0)
prob = ODEProblem(sir!, u0, tspan, p)
sol = solve(prob)

# Visualize the solution
plot(sol, label=["S" "I" "R"], legend=:right)

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

In this example, we will use [LabelledArrays.jl](https://github.com/SciML/LabelledArrays.jl) to get DSL-like syntax.
===#

using LabelledArrays
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---
function lorenz!(du,u,p,t)
    du.x = p.σ*(u.y-u.x)
    du.y = u.x*(p.ρ-u.z) - u.y
    du.z = u.x*u.y - p.β*u.z
end

#---
u0 = LVector(x=1.0, y=0.0, z=0.0)
p = LVector(σ=10.0, ρ=28.0, β=8/3)
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
sol = solve(prob)

# x-y-z time-series
plot(sol)

# `idxs=(1, 2, 3)` makes a phase plot with 1st, 2nd, and the 3rd state variable. With `LabelledArrays`, you can use symbols instead of index numbers.

plot(sol, idxs=(:x, :y, :z))

# The zeroth variable in `idxs` is the independent variable (usually time). The below command plots the time series of the second state variable (`y`).

plot(sol, idxs=(0, 2))

# ## Saving simulation results

using DataFrames
using CSV

df = DataFrame(sol)
CSV.write("lorenz.csv", df)

#===

## Using ModelingToolkit.jl (advanced)

[ModelingToolkit.jl](https://mtk.sciml.ai/dev/) is a high-level package for symbolic-numeric modeling and simulation ni the Julia DiffEq ecosystem.

===#

using DifferentialEquations
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

# ### Exponential decay model

@parameters λ       ## Decaying rate constant
@variables t C(t)   ## Time and concentration
D = Differential(t) ## Differential operator

# Define ODE equation(s)
eqs = [D(C) ~ -λ*C]

# Define an ODE system
@named expdecaySys = ODESystem(eqs)

#---
u0 = [C => 1.0]
p = [λ => 1.0]
tspan = (0.0, 2.0)

prob = ODEProblem(expdecaySys, u0, tspan, p)
sol = solve(prob)

plot(sol)

# ### SIR model

@parameters β γ
@variables t s(t) i(t) r(t)
D = Differential(t)

eqs = [
    D(s) ~ -β * s * i,
    D(i) ~ β * s * i - γ * i,
    D(r) ~ γ * i
]

@named sirSys = ODESystem(eqs)

p = [β => 1.0, γ => 0.3]
u0 = [s => 0.99, i => 0.01, r => 0.00]
tspan = (0.0, 20.0)

prob = ODEProblem(sirSys, u0, tspan, p)
sol = solve(prob)

plot(sol)

#===
## Using Catalyst.jl for chemical reaction networks

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) is a domain-specific language (DSL) package to simulate chemical reaction networks.
===#

using Catalyst
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

# ### Exponential decay model

decay_rn = @reaction_network begin
    λ, C --> 0
end

p = [:λ => 1.]
u0 = [:C => 1.]
tspan = (0., 2.)

prob = ODEProblem(decay_rn, u0, tspan, p)
sol = solve(prob)

plot(sol, title="Exponential Decay")

# ### SIR model

sir_rn = @reaction_network begin
    β, S + I --> 2I
    γ, I --> R
end

p = [:β => 1.0, :γ => 0.3]
u0 = [:S => 0.99, :I => 0.01, :R => 0.00]
tspan = (0., 20.)

prob = ODEProblem(sir_rn, u0, tspan, p)
sol = solve(prob)

plot(sol, legend=:right, title = "SIR Model")

# ## Runtime information

import Pkg
Pkg.status()

#---
import InteractiveUtils
InteractiveUtils.versioninfo()
