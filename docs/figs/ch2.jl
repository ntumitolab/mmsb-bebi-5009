# Chapter 2

# ## Fig 2.04
# Exponential decay
using Plots

pl = plot(xlabel="Time", ylabel="Concentration", title="Fig 2.4")
for k in 1:3
    plot!(pl, t -> 3 * exp(-k * t), 0, 5, label="exp(-$(k)t)")
end
pl

# ## Fig 2.09
# Metabolic network simulation
using OrdinaryDiffEq
using Catalyst

@time "Build system" rn209 = @reaction_network model209 begin
    k1, 0 --> A
    k2, A --> B
    k3, A + B --> C + D
    k4, C --> 0
    k5, D --> 0
end

pmap = (:k1 => 3.0, :k2 => 2.0, :k3 => 2.5, :k4 => 3.0, :k5 => 4.0)
u0map = (:A => 0.0, :B => 0.0, :C => 0.0, :D => 0.0)
tend = 10.0
@time "Build problem" prob209 = ODEProblem(rn209, u0map, tend, pmap)
@time "Solve problem" sol = solve(prob209, Tsit5())

# Visualize the results
plot(sol, title="Fig 2.9", xlabel="Time", ylabel="Concentration")

# ## Figure 2.11
# Intro to model reduction of ODE metabolic networks.
using OrdinaryDiffEq
using Catalyst

@time "Build system" rn211 = @reaction_network model211 begin
    k0, 0 --> A
    (k1, km1), A <--> B
    k2, B --> 0
end

pmap211 = (:k0=>0.0, :k1=>9.0, :km1=>12.0, :k2=>2.0)
u0map211 = (:A=>0.0, :B=>10.0)
tend = 3.0
@time "Build problem" prob211 = ODEProblem(rn211, u0map211, tend, pmap211)
@time "Solve problem" sol211 = solve(prob211, Tsit5())

# Fig 2.11
plot(sol211, title="Fig 2.11", xlabel="Time", ylabel="Concentration")

# ## Figure 2.12
# Rapid equilibrium assumption
@time "Build system" rn212 = @reaction_network model212 begin
    @parameters k0 k1 km1 k2
    @equations begin
        A ~ C * km1 / (km1 + k1)
        B ~ C * k1 / (km1 + k1)
    end
    k0, 0 --> C
    k2 * B, C => 0
end

tend = 3.0
u0map212 = (:C=>10.0)
pmap212 = pmap211
@time "Build problem" prob212 = ODEProblem(rn212, u0map212, tend, pmap212; mtkcompile=true)
@time "Solve problem" sol212 = solve(prob212, Tsit5())

# Fig 2.12
plot(sol211, title="Fig 2.12", xlabel="Time", ylabel="Concentration", linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(sol212, idxs=[:A, :B],linestyle=:solid, label=["A (Reduced model)" "B (Reduced model)"])

# ## Figure 2.13
# When another set of parameters is not suitable for rapid equilibrium assumption.
pmap213 = Dict(:k0=>9.0, :k1=>21.0, :km1=>12.0, :k2=>2.0)
u0map213 = (:A=>8.0, :B=>4.0)
tend = 3.0
prob213full = remake(prob211; p=pmap213, u0=u0map213, tspan=(0.0, tend))
prob213re = ODEProblem(rn212, (:C=>12.0), tend, pmap213; mtkcompile=true)

@time sol213full = solve(prob213full, Tsit5())
@time sol213re = solve(prob213re, Tsit5())

# Fig 2.13
plot(sol213full, title="Fig 2.13", xlabel="Time", ylabel="Concentration", linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(sol213re, idxs=[:A, :B], linestyle=:solid, label=["A (Reduced model)" "B (Reduced model)"])

# ## Figure 2.14 : QSSA
# Quasi-steady state assumption on species A (D(A) ~ 0)
@time "Build system" rn214 = @reaction_network model214 begin
    @parameters k0 k1 km1 k2
    @equations begin
        A ~ (k0 + km1 * B) / k1
    end
    (k1 * A, (k2 + km1)), 0 <--> B
end

pmap214 = pmap213
u0map214 = (:B=>(pmap214[:k1] * 12.0 - pmap214[:k0]) / (pmap214[:k1] + pmap214[:km1]))
@time sol214 = solve(ODEProblem(rn214, u0map214, tend, pmap214; mtkcompile=true), Tsit5())

# Fig 2.14
plot(sol213full, title="Fig 2.14", xlabel="Time", ylabel="Concentration", linestyle=:dash, label=["A (Full model)" "B (Full model)"])
plot!(sol214, idxs=[:A, :B], linestyle=:solid, label=["A (QSSA)" "B (QSSA)"])

# ## Problem 2.4.6
using OrdinaryDiffEq
using Plots

sol246 = ODEProblem((u, p, t) -> p * (1.0 - u), 0.0, 10.0, 1.0) |> solve
plot(sol246, title="Problem 2.4.6", xlabel="Time", ylabel="Concentration")
