# Chapter 1

*Cells are complex systems*: many components, interactions, and feedbacks.

## Systems

- A collection of interacting *components*.
- A *boundary* separating the system and its external environment.
- Complex Systems: **nonlinear** interactions and **feedback loops**.

### Negative feedback

- System components inhibit their own activity. e.g. a thermostat corrects for deviation of temperature from a set-point.
- Negative feedback generally stabilizes system behavior manifesting self-regulation and homeostasis.
- Instabilities and oscillations can arise when there is a lag in the action of a negative feedback loop.

### Positive feedback

- System components enhance their own activity. e.g. a microphone and an amplifier.
- Divergent (unconstrained) / saturation (constrained) behavior. e.g. Flip-flop switches.

## Models

- Models are abstractions of reality. They could be
  - physical / biological: mouse model for drug testing
  - conceptual: ball-and-stick molecular model
  - mechanistic: physico-chemical laws
  - predictive: inferred by data
- In this class, we investigate biological processes as dynamic systems with mathematical models that mimic the behavior of intracellular networks. This approach is called **Dynamic Modelling**.

## Dynamic Modelling

Focusing on **rates of change** of each component in the network.

### Aiding biological investigation

- Potential insights to mechanisms / behaviors.
- Transparent description
- Not affected by limits of lab conditions
- Omniscient and omnipotent about this system.

### To investigate mechanistic models

Simulations indicate how a system behaves; model analysis reveals why a system behaves as it does.

Model simulation:
- Direct approach: solving the *forward problem*
- Is used to predict system behavior under given conditions
- Sometimes refer to in silico experiments because they use computer to mimic the
behavior of biological systems

Model Analysis (Ch.4)
- Solving the *reverse problem* sometimes involves sophisticated mathematical techniques
- Investigate directly, yielding insight into their potential behaviors
- reveal non-intuitive connections between the structure of a system and its consequent
behavior

## BASIC FEATURES OF DYNAMIC MATHEMATICAL MODELS

### State variables

- The primary components in the system.
- The collection of all of these state variables is called the **state** of the system, a complete description of the system’s condition at any given time.
- The model’s dynamic behaviour is the time-course for the collection of state variables.

### Model parameters

- The environmental effects and interactions among system components. e.g. chemical kinetic constants, temperature, physical constants, buffered concentrations.
- Typically constant throughout the simulations.
- How to obtain
  - *Databases* like [BRENDA](https://www.brenda-enzymes.org/).
  - *Fitting* the model to the experiment results.

### Steady-State Behavior and Transient Behavior

- Steady state: persistent operating state, asymptotic behavior after a significant stretches of time.
- Transient bahavior: the immediate response of a system to perturbations.

### Linearity vs nonlinearity

- linear: $x = k_1y + k_2z$, easier to model (e.g. traditional control system models)
- nonlinear: $x = k_1y^2 + k_2z / (k_2z + 1)$. As in most biochemical models.

### Saturation

- Hyperbolic saturation: $y = \frac{x}{x + 1}$ [Graph link](https://www.desmos.com/calculator/c1nalotl5q).
- Sigmoidal saturation: $y = \frac{x^4}{x^4 + 1}$ [Graph link](https://www.desmos.com/calculator/csf4wcaqej)

### Global vs local behavior

- Over small domains, nonlinearities can always be approximated by linear relationships.
- Local approximation allows one to apply linear analysis tools. (e.g. Jacobian matrix)
- Local analysis at these points can then provide comprehensive insight into complex global behaviors.

### Deterministic vs stochastic

- Deterministic: Behavior dependent solely on model states, exactly reproducible for repeated simulations.
- Stochastic: Randomness introduced into model. Repeated simulations yield a set of results. e.g. for gene regulation networks subject to few molecule conditions sensitive to thermal noise.

## Tools to investigate dynamic mathematical models

- Differential Equations solvers: dynamic behavior
- Sensitivity Analysis: dependence of steady-state behaviour on internal and external conditions.
- Stability Analysis: phase plane analysis, characterizing long-term behaviour. (steady-state or not)
- Bifurcation Analysis: dependence of system dynamics on internal and external conditions.

## What we want to know

Intracellular dynamics
- Metabolism: chemical reaction networks, enzyme-catalysed reactions, allosteric regulation (Chapter 5)
- Signal Transduction: G protein signalling, MAPK signalling cascade, bacterial chemotaxis, calcium oscillations (Chapter 6)
- Genetic Networks: switches (lac operon, phage lambda lysis/lysogeny switch, engineered toggle switch), oscillators (Goodwin oscillator, circadian rhythms, cell cycle, repressilator), computation (Chapter 7)
- Electrophysiology: voltage-gated ion channels, Nernst potential, Morris-Lecar model, intercellular communication (gap junctions, synaptic transmission, neuronal circuits) (Chapter 8)
