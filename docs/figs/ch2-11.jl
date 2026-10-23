# # Fig 2.11-14
# Intro to model reduction of ODE metabolic networks.
using OrdinaryDiffEq
using Catalyst
using Plots

# ## Fig 2.11
@time "Build system" rn211 = @reaction_network model211 begin
    k0, 0 --> A
    (k1, km1), A <--> B
    k2, B --> 0
end

up211 = (:k0 => 0.0, :k1 => 9.0, :km1 => 12.0, :k2 => 2.0, :A => 0.0, :B => 10.0)
tend = 3.0
@time "Build problem" prob211 = ODEProblem(rn211, up211, tend)
@time "Solve problem" sol211 = solve(prob211, Tsit5())
plot(sol211, title="Fig 2.11", xlabel="Time", ylabel="Concentration")

# ## Fig 2.12
# Rapid equilibrium assumption
@time "Build system" rn212 = @reaction_network model212 begin
    @parameters k0 k1 km1 k2
    @equations begin
        A ~ C * km1 / (km1 + k1)
        B ~ C * k1 / (km1 + k1)
    end
    k0, 0 --> C
    k2 * k1 / (km1 + k1), C --> 0
end

#---
tend = 3.0
up212 = (:k0 => 0.0, :k1 => 9.0, :km1 => 12.0, :k2 => 2.0, :C => 10.0)
@time "Build problem" prob212 = ODEProblem(rn212, up212, tend; mtkcompile=true)
@time "Solve problem" sol212 = solve(prob212, Tsit5())

plot(sol211, linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(sol212, idxs=[:A, :B], linestyle=:solid, label=["A (Reduced model)" "B (Reduced model)"])
plot!(title="Fig 2.12", xlabel="Time", ylabel="Concentration")

# ## Fig 2.13
# When another set of parameters breaks rapid equilibrium assumption.
pmap213 = Dict(:k0=>9.0, :k1=>21.0, :km1=>12.0, :k2=>2.0)
u0map213 = Dict(:A=>8.0, :B=>4.0)
tend = 3.0
prob213full = remake(prob211; p=pmap213, u0=u0map213, tspan=(0.0, tend))
prob213re = ODEProblem(rn212, (:C=>12.0), tend, pmap213; mtkcompile=true)
@time sol213full = solve(prob213full, Tsit5())
@time sol213re = solve(prob213re, Tsit5())

plot(sol213full,  linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(sol213re, idxs=[:A, :B], linestyle=:solid, label=["A (Reduced model)" "B (Reduced model)"])
plot!(title="Fig 2.13", xlabel="Time", ylabel="Concentration")

# ## Fig 2.14
# Quasi-steady state assumption on species A (`D(A) ~ 0`). You can model this in `Catalyst.jl` DSL like:
@time "Build system" rn214 = @reaction_network model214 begin
    @parameters k0=9 k1=21 km1=12 k2=2
    @variables A(t)
    @equations begin
        A ~ (k0 + km1 * B) / k1
    end
    (k1 * A, (k2 + km1)), 0 <--> B
end

# Or use `ModelingToolkit.jl` equations to set `dA ~ 0` directly. The QSS concentrations of A will be calculated by `ModelingToolkit.jl`.
using ModelingToolkit

function model214(u0; name=:model214)
    @independent_variables t
    @parameters k0=9 k1=21 km1=12 k2=2

    D = Differential(t)
    @variables A(t) B(t) = (k1 * u0 - k0) / (km1 + k1)

    v0 = k0
    v1 = k1 * A - km1 * B
    v2 = k2 * B
    dA = v0 - v1
    dB = v1 - v2

    eqs = [
        0 ~ dA, ## Quasi-steady state assumption on species A (`D(A) ~ 0`)
        D(B) ~ dB,
    ]
    return ODESystem(eqs, t; name)
end

tend = 3.0
@time "Build system" @mtkcompile sys214 = model214(12.0)
@time "Build problem" prob214 = ODEProblem(sys214, [], (0.0, tend))

#---
equations(sys214)

#---
observed(sys214)
