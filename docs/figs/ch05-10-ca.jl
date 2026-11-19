# # Fig 5.10, 5.11
# Methionine model
using ComponentArrays
using OrdinaryDiffEq
using Plots
using Startup: meshgrid, normalize_gradient
Plots.gr(linewidth=1.5)

#---
hil(x, k=one(x)) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

function met_ode!(du, u, p, t)
    Hcy = u.AdoHcy * p.K_AHC / p.Adenosine                                ## Hcy ~ AdoHcy * K_AHC / Adenosine
    v_MATI = p.v_MATI_max * hil(p.Met * hil(p.K_MATI_i, u.AdoMet), p.K_MATI_m)
    K_MATIII_m1 = 20000 / (1 + 5.7 * hil(u.AdoMet, 600)^2)                ## needed by v_MATIII
    v_MATIII = p.v_MATIII_max * hil(p.Met, K_MATIII_m1 * hil(p.K_MATIII_m2, p.Met))
    v_GNMT = p.v_GNMT_max * hil(u.AdoMet, p.K_GNMT_m, 2.3) * hil(p.K_GNMT_i, u.AdoHcy)
    K_MET_m1 = 10 + 2.5 * u.AdoHcy                                        ## needed by v_MET
    v_MET = p.v_MET_max * hil(u.AdoMet, K_MET_m1) * hil(p.A_over_K_MET_m2)
    v_D = p.alpha_d * Hcy

    ## The two `D(...) ~ ...` equations of the MTK model:
    du.AdoMet = (v_MATI + v_MATIII) - (v_GNMT + v_MET)
    du.AdoHcy = (v_GNMT + v_MET - v_D) * hil(p.Adenosine, p.K_AHC)
    return nothing
end

#---
function get_gradient(f!, p, xx, yy; t = 0.0)
    dx = map(xx, yy) do x, y
        u = ComponentArray(AdoMet = x, AdoHcy = y)
        du = similar(u)
        f!(du, u, p, t)
        du.AdoMet
    end
    dy = map(xx, yy) do x, y
        u = ComponentArray(AdoMet = x, AdoHcy = y)
        du = similar(u)
        f!(du, u, p, t)
        du.AdoHcy
    end
    return dx, dy
end

# ## Figure 5.10
p510 = ComponentArray(
    K_AHC = 0.1,
    Adenosine = 1.0,
    v_MATI_max = 561.0,
    Met = 48.5,
    K_MATI_m = 41.0,
    K_MATI_i = 50.0,
    v_MATIII_max = 22870.0,
    K_MATIII_m2 = 21.1,
    v_GNMT_max = 10600.0,
    K_GNMT_m = 4500.0,
    K_GNMT_i = 20.0,
    v_MET_max = 4544.0,
    A_over_K_MET_m2 = 0.1,
    alpha_d = 1333.0,
)
u510 = ComponentArray(AdoMet = 10.0, AdoHcy = 10.0)
tend = 5.0
alg = FBDF()
@time "Build problem" prob510 = ODEProblem(met_ode!, u510, (0.0, tend), p510)
@time "Solve problem" sol510 = solve(prob510, alg)

## The MTK version plotted the whole solution, which contained exactly the
## two states `AdoMet` and `AdoHcy`.
plot(sol510, idxs=[:AdoMet, :AdoHcy], labels=["AdoMet" "AdoHcy"], title="Fig. 5.10", xlabel="Time (s)", ylabel="Concentration (μM)")

# ## Figure 5.11 A
xrange = range(0, 1200, 101)
yrange = range(0, 6, 101)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(met_ode!, prob510.p, xx, yy)

xx_sp = xx[1:5:end, 1:5:end]
yy_sp = yy[1:5:end, 1:5:end]
dx_sp = dx[1:5:end, 1:5:end]
dy_sp = dy[1:5:end, 1:5:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange) * 4 , step(yrange) * 4)
quiver(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="AdoMet nullcline")
contour!(xrange, yrange, dy, levels=[0],  line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="AdoHcy nullcline")
plot!(title="Fig. 5.11 (A)", xlabel="AdoMet (μM)", ylabel="AdoHcy (μM)", legend=:bottomright, size=(800, 600), xlims=(0, 1200), ylims=(0, 6))

# ## Figure 5.11 B
# Increased methionine level
prob511b = remake(prob510, p=[:Met => 51.0])

xrange = range(0, 1200, 101)
yrange = range(0, 6, 101)
xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(met_ode!, prob511b.p, xx, yy)

xx_sp = xx[1:5:end, 1:5:end]
yy_sp = yy[1:5:end, 1:5:end]
dx_sp = dx[1:5:end, 1:5:end]
dy_sp = dy[1:5:end, 1:5:end]
dx_sp_norm, dy_sp_norm = normalize_gradient(dx_sp, dy_sp, step(xrange) * 4 , step(yrange) * 4)

quiver(xx_sp, yy_sp, quiver=(dx_sp_norm, dy_sp_norm), color=:gray)
contour!(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [],  line=(:black, :solid), label="AdoMet nullcline")
contour!(xrange, yrange, dy, levels=[0],  line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="AdoHcy nullcline")
plot!(title="Fig. 5.11 (B)", xlabel="AdoMet (μM)", ylabel="AdoHcy (μM)", legend=:bottomright, size=(800, 600), xlims=(0, 1200), ylims=(0, 6))
