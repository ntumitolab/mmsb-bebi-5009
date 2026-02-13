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
# Get gradient over an xy grid
function get_gradient(prob, xsym, ysym; t = nothing, xrange=range(0, 2, 21), yrange=range(0, 2, 21))
    swap_or_not(x, y; xidx=1) = xidx == 1 ? [x, y] : [y, x]
    ∂F = prob.f
    ps = prob.p
    sys = prob.f.sys
    xidx = ModelingToolkit.variable_index(sys, xsym)
    yidx = ModelingToolkit.variable_index(sys, ysym)
    xx = [x for y in yrange, x in xrange]
    yy = [y for y in yrange, x in xrange]
    dx = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[xidx], xx, yy)
    dy = map((x, y) -> ∂F(swap_or_not(x, y; xidx), ps, t)[yidx], xx, yy)
    return (; xx, yy, dx, dy)
end

xrange = range(0, 2, 21)
yrange = range(0, 2, 21)
@unpack xx, yy, dx, dy = get_gradient(prob402a, s1, s2; t = nothing, xrange, yrange)

## Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.1
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength

pl44a = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), color=:gray)

for sol in sols
    plot!(pl44a, sol, idxs=(s1, s2), label=nothing)
end

plot!(pl44a, title="Fig. 4.4 A", xlabel="[S1]", ylabel="[S2]")


# ## Figure 4.5A
# Phase plot with nullcline
pl45a = contour(xrange, yrange, dx, levels=[0], cbar=false, aspect_ratio=1, size=(600, 600), color=:black)
contour!(pl45a, xrange, yrange, dy, levels=[0], cbar=false, line=(:dash, :black))
plot!(pl45a, [], [], line=(:black, :solid), label="A nullcline")
plot!(pl45a, [], [], line=(:black, :dash), label="B nullcline")
for sol in sols
    plot!(pl45a, sol, idxs=(s1, s2), label=nothing)
end
plot!(pl45a, title="Fig. 4.5 A\nPhase plot with nullclines", xlabel="[S1]", ylabel="[S2]", xlims=(0, 2), ylims=(0, 2), legend=:bottomright)


# ## Figure 4.5 B
# Vector field with nullclines
pl45b = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 2), ylims=(0, 2), color=:gray)
contour!(pl45b, xrange, yrange, dx, levels=[0], cbar=false, color=:black)
contour!(pl45b, xrange, yrange, dy, levels=[0], cbar=false, line=(:dash, :black))
plot!(pl45b, [], [], line=(:black, :solid), label="A nullcline")
plot!(pl45b, [], [], line=(:black, :dash), label="B nullcline")
