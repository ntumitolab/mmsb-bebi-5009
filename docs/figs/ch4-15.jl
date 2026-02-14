# # Fig 4.15, 4.16, 4.17
# Oscillatory networks.
# ## Figure 4.15 (A)
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)
#---
function model415(; tend=8.0)
    @parameters k0=8.0 k1=1.0 k2=5.0 n=2.0
    @variables A(t)=1.5 B(t)=1.0
    vAB = k1 * A * (1 + B^n)
    eqs = [
        D(A) ~ k0 - vAB
        D(B) ~ vAB - k2 * B
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end
@time prob415 = model415()

#---
@unpack A, B = prob415.f.sys
u0s = [
    [A=>1.5, B=>1.0],
    [A=>0.0, B=>1.0],
    [A=>0.0, B=>3.0],
    [A=>2.0, B=>0.0],
]

@time sols = map(u0s) do u0
    solve(remake(prob415, u0=u0), Tsit5())
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.15 A")

# ## Fig 4.15 (B)
# Vector field with nullclines
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

@unpack xx, yy, dx, dy = get_gradient(prob415, A, B; t = nothing, xrange=range(0, 4, 21), yrange=range(0, 4, 21))

## Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.7
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength

pl415b = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 4), ylims=(0, 4), color=:gray)

## Finer resolution for nullclines
xrange = range(0, 4, 51)
yrange = range(0, 4, 51)
@unpack xx, yy, dx, dy = get_gradient(prob415, A, B; t = nothing, xrange, yrange)

contour!(pl415b, xrange, yrange, dx, levels=[0], cbar=false, line=(:black))
contour!(pl415b, xrange, yrange, dy, levels=[0], cbar=false, line=(:dash, :black))
plot!(pl415b, [], [], line=(:black, :solid), label="A nullcline")
plot!(pl415b, [], [], line=(:black, :dash), label="B nullcline")
plot!(pl415b, title="Fig 4.15 B", xlabel="[A]", ylabel="[B]", legend=:topright)
pl415b

#===
## Fig 4.16 A

Oscillatory parameter set
===#
@unpack n = prob415.f.sys
prob416 = remake(prob415, p=[n=>2.5])
u0s = [
    [A=>1.5, B=>1.0],
    [A=>0.0, B=>1.0],
    [A=>0.0, B=>3.0],
    [A=>2.0, B=>0.0],
]

@time sols = map(u0s) do u0
    solve(remake(prob416, u0=u0), Tsit5())
end

plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.16 A")

# ## Fig 4.16 b
@unpack xx, yy, dx, dy = get_gradient(prob416, A, B; t = nothing, xrange=range(0, 4, 21), yrange=range(0, 4, 21))

## Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.7
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength
pl416b = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 4), ylims=(0, 4), color=:gray)

fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.16 B\nPhase plot",
    aspect = 1,
)

streamplot!(ax, ∂F416, xx, yy)
contour!(ax, xx, yy, ∂A416, levels=[0], color=:black, linestyle=:solid, linewidth=2, label="A nullcline")
contour!(ax, xx, yy, ∂B416, levels=[0], color=:black, linestyle=:dash, linewidth=2, label="B nullcline")
for sol in sols
    aa = [sol(t).A for t in ts]
    bb = [sol(t).B for t in ts]
    lines!(ax, aa, bb, color=:tomato)
end
axislegend(ax, position = :rc)
limits!(ax, 0.0, 4.0, 0.0, 4.0)
fig

# ## Fig 4.17
prob417 = remake(prob415, p=ps416, u0=ComponentArray(A=2.0, B=1.5), tspan=(0.0, 10.0))
@time sol = solve(prob417, Tsit5())

fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    xlabel = "[A]",
    ylabel = "[B]",
    title = "Fig. 4.17",
    aspect = 1,
)

lines!(ax, sol, idxs=(1, 2), color=:tomato, label="Trajectory")
xx = 1:0.01:3
yy = 1:0.01:3
∂A417 = [_dA415(x, y, ps416, nothing) for x in xx, y in yy]
∂B417 = [_dB415(x, y, ps416, nothing) for x in xx, y in yy];

streamplot!(ax, ∂F416, xx, yy)
contour!(ax, xx, yy, ∂A417, levels=[0], color=:black, linestyle=:solid, linewidth=2, label="A nullcline")
contour!(ax, xx, yy, ∂B417, levels=[0], color=:black, linestyle=:dash, linewidth=2, label="B nullcline")
axislegend(ax, position = :rc)
limits!(ax, 1.0, 3.0, 1.0, 3.0)
fig
