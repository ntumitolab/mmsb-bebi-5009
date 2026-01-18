# # Fig 7.21
# repressilator model
using OrdinaryDiffEq
using ComponentArrays: ComponentArray
using SimpleUnPack
using CairoMakie

#---
function model721!(D, u, p, t)
    @unpack α, α0, n, β = p
    @unpack mA, mB, mC, A, B, C = u
    D.mA = α / (1 + C^n) + α0 - mA
    D.mB = α / (1 + A^n) + α0 - mB
    D.mC = α / (1 + B^n) + α0 - mC
    D.A = β * (mA - A)
    D.B = β * (mB - B)
    D.C = β * (mC - C)
    nothing
end

#---
ps721 = ComponentArray(
    α = 298.2,
    α0 = 0.03,
    n = 2.0,
    β = 0.2
)
u0721 = ComponentArray(
    mA = 0.2,
    mB = 0.3,
    mC = 0.4,
    A = 0.1,
    B = 0.1,
    C = 0.5
)

tspan = (0.0, 400.0)

prob721 = ODEProblem(model721!, u0721, tspan, ps721)

#---
@time sol721 = solve(prob721, KenCarp47())
#---
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time",
    ylabel = "Concentration",
    title = "Fig 7.21"
)
lines!(ax, 0..400, t-> sol721(t).A, label = "A")
lines!(ax, 0..400, t-> sol721(t).B, label = "B")
lines!(ax, 0..400, t-> sol721(t).C, label = "C")
axislegend(ax, position = :lt)
fig
