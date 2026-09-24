# # Fig 7.21
# Repressilator model.
using DifferentialEquations
using Catalyst
using Plots
Plots.gr(linewidth=1.5)

@time "Build system" rn721 = @reaction_network begin
    (α0 + hillr(C, α, 1, n), 1), 0 <--> mA
    (α0 + hillr(A, α, 1, n), 1), 0 <--> mB
    (α0 + hillr(B, α, 1, n), 1), 0 <--> mC
    (β * mA, β), 0 <--> A
    (β * mB, β), 0 <--> B
    (β * mC, β), 0 <--> C
end

#---
alg = FBDF()
up721 = Dict(
    :α => 298.2, :α0 => 0.03, :n => 2.0, :β => 0.2,
    :mA => 0.2, :mB => 0.3, :mC => 0.4,
    :A => 0.1, :B => 0.1, :C => 0.5
)

tend = 300.0
@time "Build problem" prob721 = ODEProblem(rn721, up721, (0.0, tend))
@time "Solve problem" sol721 = solve(prob721, alg)

plot(sol721, idxs=[:A, :B, :C], xlabel="Time", ylabel="Concentration", title="Fig 7.21", legend=:left)
