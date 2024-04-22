#===
# Fig 1.07

Collins toggle switch

For Figures 1.7, 7.13, 7.14, 7.15
===#

using DifferentialEquations
using ModelingToolkit
using Plots
Plots.default(linewidth=2)

# Convenience functions
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)
exprel(x) = x / expm1(x)

# Define the problem
function build_collins(;name)
    @parameters begin
        a1=3.0
        a2=2.5
        β=4.0
        γ=4.0
    end

    @variables begin
        t
        s1(t)=0.075
        s2(t)=2.5
        i1(t)
        i2(t)
    end

    D = Differential(t)

    eqs = [
        D(s1) ~ a1 * hil(1 + i2, s2, β) - s1
        D(s2) ~ a2 * hil(1 + i1, s1, γ) - s2
        i2 ~ 10 * (10 < t) * (t < 20)
        i1 ~ 10 * (30 < t) * (t < 40)
    ]
    sys = ODESystem(eqs; name)
    return structural_simplify(sys)
end

#---
@named sys = build_collins()

# Solve the problem
tspan = (0.0, 50.0)
prob = ODEProblem(sys, [], tspan, [])
sol = solve(prob, tstops=10:10:40)

# Visual
plot(sol, legend=:right, xlabel = "Time", ylabel="Concentration", title="Fig 1.7")
