# # Fig 4.7, 4.8
# Symmetric (bistable) biological networks.
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)

#---
function model407(; tend = 4.0)
    @parameters k1=20 k2=20 k3=5 k4=5 n1=4 n2=1
    @variables A(t)=3 B(t)=1
    eqs = [
        D(A) ~ k1 / (1 + B^n1) - k3 * A
        D(B) ~ k2 / (1 + A^n2) - k4 * B
    ]
    @mtkcompile sys = ODESystem(eqs, t)
    return ODEProblem(sys, [], tend)
end

#===
## Fig 4.7 A

Asymmetric parameter set
===#
@time prob407 = model407()

#---
@unpack A, B = prob407.f.sys
@time sol1 = solve(prob407, Tsit5())
@time sol2 = solve(remake(prob407, u0=[A=>1.0, B=>3.0]), Tsit5())
pl47_1 = plot(sol1, xlabel="Time", ylabel="Concentration", title= "Fig 4.7A")
pl47_2 = plot(sol2, xlabel="Time", ylabel="Concentration")
plot(pl47_1, pl47_2, layout=(2, 1))

# ## Fig 4.7 B
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

@unpack xx, yy, dx, dy = get_gradient(prob407, A, B; t = nothing, xrange=range(0, 5, 21), yrange=range(0, 5, 21))
## Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.5
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength

pl47b = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 5), ylims=(0, 5), color=:gray)

## Finer resolution for nullclines
xrange = range(0, 5, 51)
yrange = range(0, 5, 51)
@unpack xx, yy, dx, dy = get_gradient(prob407, A, B; t = nothing, xrange, yrange)
contour!(pl47b, xrange, yrange, dx, levels=[0], cbar=false, line=(:black))
contour!(pl47b, xrange, yrange, dy, levels=[0], cbar=false, line=(:dash, :black))
plot!(pl47b, [], [], line=(:black, :solid), label="A nullcline")
plot!(pl47b, [], [], line=(:black, :dash), label="B nullcline")
plot!(pl47b, title="Fig 4.7 B", xlabel="[A]", ylabel="[B]", legend=:topright)
pl47b

# ## Fig 4.8
# Symmetric parameter set
@unpack k1, k2, k3, k4, n1, n2 = prob407.f.sys
prob408 = remake(prob407, p=[k1=>20, k2=>20, k3=>5, k4=>5, n1=>4, n2=>4])

@time sol1 = solve(prob408, Tsit5())
@time sol2 = solve(remake(prob408, u0=[A=>1.0, B=>3.0]), Tsit5())

pl48_1 = plot(sol1, xlabel="Time", ylabel="Concentration", title= "Fig. 4.8 A")
pl48_2 = plot(sol2, xlabel="Time", ylabel="Concentration")
plot(pl48_1, pl48_2, layout=(2, 1))

#---
@unpack xx, yy, dx, dy = get_gradient(prob408, A, B; t = nothing, xrange=range(0, 5, 21), yrange=range(0, 5, 21))
## Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.5
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength

pl48b = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 5), ylims=(0, 5), color=:gray)

## Finer resolution for nullclines
xrange = range(0, 5, 51)
yrange = range(0, 5, 51)
@unpack xx, yy, dx, dy = get_gradient(prob408, A, B; t = nothing, xrange, yrange)
contour!(pl48b, xrange, yrange, dx, levels=[0], cbar=false, line=(:black))
contour!(pl48b, xrange, yrange, dy, levels=[0], cbar=false, line=(:dash, :black))
plot!(pl48b, [], [], line=(:black, :solid), label="A nullcline")
plot!(pl48b, [], [], line=(:black, :dash), label="B nullcline")
plot!(pl48b, title="Fig 4.8 B", xlabel="[A]", ylabel="[B]", legend=:topright)
pl48b

# ## Fig 4.8 C
# Around the unstable steady-state
@unpack xx, yy, dx, dy = get_gradient(prob408, A, B; t = nothing, xrange=range(1.0, 1.5, 21), yrange=range(1.0, 1.5, 21))
## Normalize vector field
maxnorm = maximum(hypot.(dx, dy))
maxlength=0.05
dxnorm = @. dx / maxnorm * maxlength
dynorm = @. dy / maxnorm * maxlength
pl48c = quiver(xx, yy, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(1, 1.5), ylims=(1, 1.5), color=:gray)

contour!(pl48c, range(1.0, 1.5, 21), range(1.0, 1.5, 21), dx, levels=[0], cbar=false, line=(:black))
contour!(pl48c, range(1.0, 1.5, 21), range(1.0, 1.5, 21), dy, levels=[0], cbar=false, line=(:dash, :black))
plot!(pl48c, [], [], line=(:black, :solid), label="A nullcline")
plot!(pl48c, [], [], line=(:black, :dash), label="B nullcline")
plot!(pl48c, title="Fig 4.8 C", xlabel="[A]", ylabel="[B]", legend=:topright)
pl48c


# ## Nullclines vs parameters
# Another way to draw nullclines is to find the analytical solutions when dA (or dB) is zero. And then sketch the nullclines in a parameteric plot.
nca47(b, p) = p.k1 / p.k3 / (1 + b^p.n1)
ncb47(a, p) = p.k2 / p.k4 / (1 + a^p.n2)

axs = map((8.0, 16.0, 20.0, 35.0)) do k1
    ps = (k1=k1, k2=20., k3=5., k4=5., n1=4., n2=4.)
    pl = plot(t->nca47(t, ps), identity, 0, 7, color=:red, label="Nullcline A")
    plot!(pl, identity, t->ncb47(t, ps), 0, 7, color=:blue, label="Nullcline B")
    plot!(pl, title="K1 = $k1", xlabel="[A]", ylabel="[B]", xlims=(0, 7), ylims=(0, 7))
    pl
end

plot(axs..., layout=(2, 2), size=(800, 800))
