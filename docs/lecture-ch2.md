# Chapter 2

**Modeling of Chemical Reaction Network**

## Chemical reaction network

The foundation for dynamic modeling of cellular behavior.

- Component (state variables): molecular species, chemically distinct molecules
- Interaction (arrows): a range of process forming a chemical reaction network.
  - Binding/unbinding
  - Reaction catalysis
  - Activity regulation
- Rates: production, consumption, and transport of the molecules

The organization of the network is apparent if we re-arrange the reactions in the form of an interaction graph.

### Closed vs open network

- Closed system: No exchange of material  between the system and the environment. At steady state, all net reaction rates are zero, i.e. thermal equilibrium.
- Open system: exchanging material with outside environment. At steady state, there are steady flux throght the system. Most biochemical networks belong to this.

### Assumptions about Chemical reaction network

So that we could use ODE systems to solve the problem.

- *Spatial homogeneity*: Well stirred. Otherwise, we should use compartmentalizations or even PDEs. (e.g. action potentials travelling along the axon)
- *Continuum hypothesis*: Great many molecules to smooth out random noise. Otherwise, we should use stochastic simulation. (e.g. gene regulatory networks)

## The Law of Mass Action

The rate of a chemical reaction is proportional to the product of the concentrations of the reactants.
- On the molecular level: the probability of a reaction occurring is proportional to the probability of the reactants colliding with one another.
- Kinetic order: the *exponent* on the reactant in a reaction.
- Rate constant

## Using derivatives to describe rates of change

### Decay reaction

For a decay reaction `A => 0`, the rate of change of concentration of species A could be described as

$$
\frac{da(t)}{dt} = -ka(t)
$$

The solution is a exponential decay: ($A_0$ is the initial concentration of A)

$$
a(t) = A_0e^{-kt}
$$

- Time constant $\tau = 1/k$ (also called the time scale, telling how fast will the system go)
- Half-life = $\tau_{0.5} = \tau ln2$

### Production and decay

For a chemical reaction `0 => A => 0`, with production rate $k_0$ and decaying rate $k_1$, the rate of change of concentration of species A could be described as

$$
\frac{da(t)}{dt} = k_0 - k_1a(t)
$$

The solution is an asymptotic to the steady-state concentration of $k_0/k_1$.

$$
a(t) = (A_0 - k_0/k_1) e^{-k_1t} + k_0/k_1
$$

when $t \rightarrow \infty$, $a(t) \rightarrow k_0/k_1$.

### Irreversible Conversion in closed system

For a reaction `A => B`, the rates of change of concentrations of species A and B could be described as

$$
\begin{aligned}
\frac{d}{dt}a(t) &= -ka(t)  \\
\frac{d}{dt}b(t) &= ka(t)  \\
\end{aligned}
$$

And the solution is

$$
\begin{aligned}
a(t) &= A_0e^{-kt}  \\
b(t) &= B_0 + A_0 - A_0e^{-kt} = B_0 + A_0 - a(t)\\
\end{aligned}
$$

The system can be reduced to a single ODE by using the *conservation relationship* $A_0 + B_0 = const$.

### Reversible Conversion in closed system

For a reaction `A <=> B`, the rates of change of concentrations of species A and B could be described as

$$
\begin{aligned}
\frac{d}{dt}a(t) &= -ka(t) + k_{-1}b(t)  \\
\frac{d}{dt}b(t) &= ka(t) - k_{-1}b(t)   \\
\end{aligned}
$$

And the solution is

$$
\begin{aligned}
a(t) &= (A_0 - a_{ss})e^{-(k_1 + k_{-1})t} + a_{ss}  \\
b(t) &= T - a(t) \\
a_{ss} &= \frac{k_{-1}}{k_{-1} + k_1}T \\
T &= A_0 + B_0
\end{aligned}
$$

The equilibrium constant $K_{eq} = \frac{b_{ss}}{a_{ss}} = \frac{k_1}{k_{-1}}$

## Differential Equations

- Ordinary differential equation (ODE): species vs time
- Partial differential equation (PDE): species vs time and space
- Autonomous differential equation: does not explicitly depend on the independent variable (e.g. time).
- Non-autonomous: some term(s) rely on the independent variable. (e.g. some parameters changes with time)

## Numerial simulations of differential equations

Proceeding step by step, treating the solution trajectory as linear. The higher order methods usually better eliminate truncation errors in the Taylor series.
- (Forward) Euler method
- Runge-Kutta methods e.g. (RK4)
- Dormand-Prince (DP5, `ode45` in MATLAB): adpative. The solver will find `dt` for itself.
- Tsitouras (`Tsit5` in Julia): a more efficient DP5.

These [solvers](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods) above use Butcher tableau to describe how to calculate the intermediate steps and finally the next step.
