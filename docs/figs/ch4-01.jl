#===
# Fig 4.2, 4.3, 4.4, 4.5

Steady states and phase plots in an asymmetric network.
===#
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)

#---
function model401(; tend = 1.5)
    @parameters k1=20 k2=5 k3=5 k4=5 k5=2 n=4
    @variables s1(t)=0.0 s2(t)=0.0
    eqs = [
        D(s1) ~ k1 / (1 + s2^n) - (k3 + k5) * s1
        D(s2) ~ k2 + k5 * s1 - k4 * s2
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

# ## Fig 4.2 A
@time prob402a = model401()

@unpack s1, s2 = prob402a.f.sys

u0s = [
    [s1 => 0.0, s2 => 0.0],
    [s1 => 0.5, s2 => 0.6],
    [s1 => 0.17, s2 => 1.1],
    [s1 => 0.25, s2 => 1.9],
    [s1 => 1.85, s2 => 1.7]
]

#---
@time sols = map(u0s) do u0
    solve(remake(prob402a, u0=u0), Tsit5())
end

#---
plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2 A")


# ## Fig. 4.2 B (Phase plot)
plot(sols[1], idxs=(s1, s2), xlabel="[S1]", ylabel="[S2]", title="Fig. 4.2 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), legend=nothing)


# ## Fig. 4.3 A (Multiple time series)
pl43 = plot(xlabel = "Time", ylabel = "Concentration", title = "Fig. 4.3 A")

for sol in sols
    plot!(pl43, sol, alpha=0.7, label=nothing)
end

pl43

# ## Fig. 4.3 B (Phase plot)
pl43B = plot(xlabel="[S1]", ylabel="[S2]", title = "Fig. 4.3 B", aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), legend=nothing)

for sol in sols
    plot!(pl43B, sol, idxs=(s1, s2), label=nothing)
end

pl43B
# ## Fig. 4.4 A
# Vector fields in phase plots
function plot_gradient(prob, xsym, ysym; maxlength=1, xmin=0, ymin=0, xmax=2, ymax=2, xstep=0.1, ystep=0.1, t = nothing, aspect_ratio=1, size=(600, 600), kwargs...)
    swap_or_not(x, y; xidx=1) = xidx == 1 ? [x, y] : [y, x]

    ∂F = prob.f
    ps = prob.p
    sys = prob.f.sys
    xidx = ModelingToolkit.variable_index(sys, xsym)
    yidx = ModelingToolkit.variable_index(sys, ysym)
    ## Raw vector field
    xx = [x for y in ymin:ystep:ymax, x in xmin:xstep:xmax]
    yy = [y for y in ymin:ystep:ymax, x in xmin:xstep:xmax]
    dx = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[xidx], xx, yy)
    dy = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[yidx], xx, yy)
    ## Normalize vector field
    maxnorm = maximum(hypot.(dx, dy))
    @. dx = dx / maxnorm * maxlength
    @. dy = dy / maxnorm * maxlength
    quiver(xx, yy, quiver=(dx, dy); aspect_ratio=aspect_ratio, size=size, xlims=(xmin, xmax), ylims=(ymin, ymax), kwargs...)
end

pl44a = plot_gradient(prob402a, s1, s2; maxlength=0.1, color=:lightblue)

for sol in sols
    plot!(pl44a, sol, idxs=(s1, s2), color=:black, label=nothing)
end

plot!(pl44a, title="Fig. 4.4 A", xlabel="[S1]", ylabel="[S2]")


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
