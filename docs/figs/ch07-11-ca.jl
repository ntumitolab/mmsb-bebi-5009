# # Fig 7.11
# Model of phage lambda decision switch.
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: meshgrid, normalize_gradient
Plots.gr(linewidth=1.5)

#---
function lambda_ode!(du, u, p, t)
    rd = u.r / 2
    cd = u.c / 2
    f1 = p.K1 * rd^2
    f2 = p.K2 * rd
    f3 = p.K3 * cd
    f4 = p.K4 * cd
    den = 1 + f1 * (1 + f2) + f3 * (1 + f4)
    du.r = p.a * (1 + 10 * f1) / den - p.delta_r * u.r
    du.c = p.b * (1 + f3) / den - p.delta_c * u.c
    return nothing
end

function get_gradient(f!, p, xname, yname, xx, yy; t = 0.0)
    du = ComponentArray(r = 0.0, c = 0.0)
    u = ComponentArray(r = 0.0, c = 0.0)
    v = map(xx, yy) do x, y
        u[xname] = x
        u[yname] = y
        f!(du, u, p, t)
        (du[xname], du[yname])
    end
    dx = map(first, v)
    dy = map(last, v)
    return (dx, dy)
end

#---
p = ComponentArray(K1=1.0, K2=0.1, K3=5.0, K4=0.5, delta_r=0.02, delta_c=0.02, a=5.0, b=50.0)
u0 = ComponentArray(r=0.0, c=0.0)

@time "Build problem" prob711 = ODEProblem(lambda_ode!, u0, (0.0, 6000.0), p)

# ## Fig 7.11 (A)
xrange = range(0, 250, 101)
yrange = range(0, 250, 101)

xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(lambda_ode!, prob711.p, :r, :c, xx, yy)
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
# Faster degradation of the cI repressor removes the lysogenic branch.
prob711b = remake(prob711; p=[:delta_r => 0.2])
dx, dy = get_gradient(lambda_ode!, prob711b.p, :r, :c, xx, yy)
dxnorm, dynorm = normalize_gradient(dx[1:5:end, 1:5:end], dy[1:5:end, 1:5:end], 5*step(xrange), 5*step(yrange))
quiver(xx_sp, yy_sp, quiver=(dxnorm, dynorm); aspect_ratio=1, size=(600, 600), xlims=(0, 250), ylims=(0, 250), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="R nullcline")
contour!(xrange, yrange, dy, levels=[0],  line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="C nullcline")
plot!(xlabel="[cI] (nM)", ylabel="[cro] (nM)", title="Fig. 7.11 (B)", legend=:topright)
