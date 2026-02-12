#===
# Fig 4.2, 4.3, 4.4, 4.5

Steady states and phase plots in an asymmetric network.
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box)
import PythonPlot as plt

#---
function model401(; tend = 1.5)
    @parameters k1=20 k2=5 k3=5 k4=5 k5=2 n=4
    @variables A(t)=0.0 B(t)=0.0
    eqs = [
        D(A) ~ k1 / (1 + B^n) - (k3 + k5) * A
        D(B) ~ k2 + k5 * A - k4 * B
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

# ## Fig 4.2 A
@time prob402a = model401()

@unpack A, B = prob402a.f.sys

u0s = [
    [A => 0.0, B => 0.0],
    [A => 0.5, B => 0.6],
    [A => 0.17, B => 1.1],
    [A => 0.25, B => 1.9],
    [A => 1.85, B => 1.7]
]

#---
@time sols = map(u0s) do u0
    solve(remake(prob402a, u0=u0), Tsit5())
end

#---
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2 A")
fig

# ## Fig. 4.2 B (Phase plot)
plot(sols[1], idxs=(A, B), xlabel="[A]", ylabel="[B]", title="Fig. 4.2 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), legend=nothing)


# ## Fig. 4.3 A (Multiple time series)
pl43 = plot(xlabel = "Time", ylabel = "Concentration", title = "Fig. 4.3 A")

for sol in sols
    plot!(pl43, sol, alpha=0.7, label=nothing)
end

pl43

# ## Fig. 4.3 B (Phase plot)
pl43B = plot(xlabel="[A]", ylabel="[B]", title = "Fig. 4.3 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), legend=nothing)

for sol in sols
    plot!(pl43B, sol, idxs=(A, B), label=nothing)
end

pl43B
# ## Fig. 4.4 A
# Vector fields in phase plots
∂F44 = prob402a.f
ps402a = prob402a.p
xx = collect(0:0.1:2)
yy = collect(0:0.1:2)
## X and Y velocities are reversed for streamplot in matplotlib
U = [∂F44([x, y], ps402a, nothing)[2] for x in xx, y in yy]
V = [∂F44([x, y], ps402a, nothing)[1] for x in xx, y in yy]

plt.figure(figsize=(6, 6))
plt.quiver(xx, yy, U, V, hypot.(U, V))
for sol in sols
    ts = 0:0.01:1.5
    aa = sol(ts, idxs=A)
    bb = sol(ts, idxs=B)
    plt.plot(aa, bb, color="black")
end
plt.xlabel("[A]")
plt.ylabel("[B]")
plt.title("Fig. 4.4 A\nVector field in phase plot")
plt.xlim(0, 2)
plt.ylim(0, 2)
plt.gcf()


# ## Figure 4.5A
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.5 A\nPhase plot with nullclines",
    aspect = 1,
)

## Phase plots
for sol in sols
    ts = 0:0.01:tend
    aa = [sol(t).A for t in ts]
    bb = [sol(t).B for t in ts]
    lines!(ax, aa, bb, color=:black)
end

## nullclines
xs = 0:0.01:2
ys = 0:0.01:2
zA44 = [_dA401(x, y, ps402a, nothing) for x in xs, y in ys]
zB44 = [_dB401(x, y, ps402a, nothing) for x in xs, y in ys]
contour!(ax, xs, ys, zA44, levels=[0], color=:black, linestyle=:solid, linewidth=2, label="A nullcline")
contour!(ax, xs, ys, zB44, levels=[0], color=:black, linestyle=:dash, linewidth=2, label="B nullcline")
limits!(ax, 0.0, 2.0, 0.0, 2.0)
axislegend(ax, position = :rb)
fig

# ## Figure 4.5 B
# Vector field
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.5 B\nVector field with nullclines",
    aspect = 1,
)
## Nullclines
contour!(ax, xs, ys, zA44, levels=[0], color=:black, linestyle=:solid, linewidth=2, label="A nullcline")
contour!(ax, xs, ys, zB44, levels=[0], color=:black, linestyle=:dash, linewidth=2, label="B nullcline")
streamplot!(ax, ∂F44, 0..2, 0..2)
limits!(ax, 0.0, 2.0, 0.0, 2.0)
fig
