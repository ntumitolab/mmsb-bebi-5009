# # Fig 5.10, 5.11
# Methionine model
using Startup: meshgrid, get_gradient, normalize_gradient
using OrdinaryDiffEq
using ModelingToolkit
using Plots
Plots.gr(linewidth=1.5)

#---
function model510(; name=:model510)
    hil(x, k=one(x)) = x / (x + k)
    hil(x, k, n) = hil(x^n, k^n)
    @independent_variables t
    D = Differential(t)
    @parameters begin
        K_AHC= 0.1
        Adenosine = 1.0
        v_MATI_max = 561.0
        Met = 48.5
        K_MATI_m = 41.0
        K_MATI_i = 50.0
        v_MATIII_max = 22870.0
        K_MATIII_m2 = 21.1
        v_GNMT_max = 10600.0
        K_GNMT_m = 4500.0
        K_GNMT_i = 20.0
        v_MET_max = 4544.0
        A_over_K_MET_m2 = 0.1
        alpha_d = 1333.0
    end
    @variables AdoMet(t)=10 AdoHcy(t)=10
    @variables Hcy(t) v_MATI(t) v_MATIII(t) K_MATIII_m1(t) v_GNMT(t) K_MET_m1(t) v_MET(t) v_D(t)

    eqs = [
        Hcy ~ AdoHcy * K_AHC / Adenosine,
        v_MATI ~ v_MATI_max * hil(Met * hil(K_MATI_i, AdoMet), K_MATI_m),
        v_MATIII ~ v_MATIII_max * hil(Met, K_MATIII_m1 * hil(K_MATIII_m2, Met)),
        K_MATIII_m1 ~ 20000 / (1 + 5.7 * hil(AdoMet, 600)^2),
        v_GNMT ~ v_GNMT_max * hil(AdoMet, K_GNMT_m, 2.3) * hil(K_GNMT_i, AdoHcy),
        K_MET_m1 ~ 10 + 2.5 * AdoHcy,
        v_MET ~ v_MET_max * hil(AdoMet, K_MET_m1) * hil(A_over_K_MET_m2),
        v_D ~ alpha_d * Hcy,
        D(AdoMet) ~ (v_MATI + v_MATIII) - (v_GNMT + v_MET),
        D(AdoHcy) ~ (v_GNMT + v_MET - v_D) * hil(Adenosine, K_AHC)
    ]

    sys = ODESystem(eqs, t; name)
    return mtkcompile(sys; simplify)
end

# ## Figure 5.10
@time "Build system" sys = model510()
tend = 5.0
alg = FBDF()
@time "Build problem" prob510 = ODEProblem(sys, [], (0.0, tend))
@time "Solve problem" sol510 = solve(prob510, alg)

plot(sol510, title="Fig. 5.10", xlabel="Time (s)", ylabel="Concentration (μM)")

# ## Figure 5.11 A
xrange = range(0, 1200, 101)
yrange = range(0, 6, 101)
xx, yy = meshgrid(xrange, yrange)
@unpack AdoMet, AdoHcy = prob510.f.sys
dx, dy = get_gradient(prob510, AdoMet, AdoHcy, xx, yy)

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
prob511b = remake(prob510, p=[sys.Met => 51.0])

xrange = range(0, 1200, 101)
yrange = range(0, 6, 101)
xx, yy = meshgrid(xrange, yrange)
@unpack AdoMet, AdoHcy = prob511b.f.sys
dx, dy = get_gradient(prob511b, AdoMet, AdoHcy, xx, yy)

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
