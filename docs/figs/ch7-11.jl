# # Fig 7.11
# Model of phage lambda decision switch
using Startup: meshgrid, get_gradient, normalize_gradient
using DifferentialEquations
using SteadyStateDiffEq
using DiffEqCallbacks
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using Plots
Plots.gr(linewidth=1.5)

function build_sys711(; name=:model711, simplify=true)
    @parameters K1=1 K2=0.1 K3=5 K4=0.5 delta_r=0.02 delta_c=0.02 a=5 b=50
    @variables r(t)=0 c(t)=0 rd(t) cd(t)
    rd = r / 2
    cd = c / 2
    f1 = K1 * rd^2
    f2 = K2 * rd
    f3 = K3 * cd
    f4 = K4 * cd
    den = 1 + f1 * (1 + f2) + f3 * (1 + f4)
    Dr = a * (1 + 10 * f1) / den - delta_r * r
    Dc = b * (1 + f3) / den - delta_c * c
    eqs = [
        D(r) ~ Dr,
        D(c) ~ Dc
    ]
    sys = ODESystem(eqs, t; name)
    return simplify ? mtkcompile(sys) : sys
end

#---
@time "Build system" sys711 = build_sys711()
@time "Build problem" prob711 = ODEProblem(sys711, [], (0.0, 6000.0))

# ## Fig 7.11 (A)
xrange = range(0, 250, 101)
yrange = range(0, 250, 101)

@unpack r, c, delta_r = sys711
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(prob711, r, c, xx, yy)
xx_sp = xx[1:5:end, 1:5:end]
yy_sp = yy[1:5:end, 1:5:end]
dxnorm, dynorm = normalize_gradient(dx[1:5:end, 1:5:end], dy[1:5:end, 1:5:end], 5*step(xrange), 5*step(yrange))

quiver(xx_sp, yy_sp, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 250), ylims=(0, 250), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="R nullcline")
contour!(xrange, yrange, dy, levels=[0],  line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="C nullcline")
plot!(xlabel="[cI] (nM)", ylabel="[cro] (nM)", title="Fig. 7.11 (A)", legend=:topright)

# ## Fig 7.11 (B)
prob711b = remake(prob711, p=[delta_r => 0.2])
dx, dy = get_gradient(prob711b, r, c, xx, yy)
dxnorm, dynorm = normalize_gradient(dx[1:5:end, 1:5:end], dy[1:5:end, 1:5:end], 5*step(xrange), 5*step(yrange))
quiver(xx_sp, yy_sp, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 250), ylims=(0, 250), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="R nullcline")
contour!(xrange, yrange, dy, levels=[0],  line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="C nullcline")
plot!(xlabel="[cI] (nM)", ylabel="[cro] (nM)", title="Fig. 7.11 (B)", legend=:topright)
