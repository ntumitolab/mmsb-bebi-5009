#===
# Fig 7.11

Model of phage lambda decision switch
===#
using OrdinaryDiffEq
using ComponentArrays
using SimpleUnPack
using Plots
Plots.default(linewidth=2)

# Model
function model711!(D, u, p, t)
    @unpack K1, K2, K3, K4, delta_r, delta_c, a, b = p
    @unpack r, c = u
    rd = r / 2
    cd = c / 2
    f1 = K1 * rd^2
    f2 = K2 * rd
    f3 = K3 * cd
    f4 = K4 * cd
    den = 1 + f1 * (1 + f2) + f3 * (1 + f4)
    D.r = a * (1 + 10 * f1) / den - delta_r * r
    D.c = b * (1 + f3) / den - delta_c * c
    return (; dr = D.r, dc = D.c)
end

#---
ps711 = ComponentArray(
    K1=1,
    K2=0.1,
    K3=5,
    K4=0.5,
    delta_r=0.02,
    delta_c=0.02,
    a=5,
    b=50,
)

u0711 = ComponentArray(
    r=0.0,
    c=0.0,
)

# ## Fig 7.11 (A)
tend = 6000.0
prob = ODEProblem(model711!, u0711, tend, ps711)

∂R = (x, y) -> model711!(ComponentArray(r=0.0, c=0.0), ComponentArray(r=x, c=y), ps711, nothing).dr
∂C = (x, y) -> model711!(ComponentArray(r=0.0, c=0.0), ComponentArray(r=x, c=y), ps711, nothing).dc
∂F = function (x, y; scale=0.2)
    dR, dC = model711!(ComponentArray(r=0.0, c=0.0), ComponentArray(r=x, c=y), ps711, nothing)
    s = sqrt(hypot(dR, dC)) * scale
    return (dR / s, dC / s)
end

rx = range(0, 250, 251)
ry = range(0, 250, 251)
rxy = range(0, 250, 21)
xx = [x for y in rxy, x in rxy]
yy = [y for y in rxy, x in rxy]

fig = plot(title="Figure 7.11 (A)")
contour!(fig, rx, ry, ∂R, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="R nullcline")
contour!(fig, rx, ry, ∂C, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="C nullcline")
quiver!(fig, xx, yy, quiver=∂F, line=(:lightblue), arrow=(:closed))

plot!(fig, xlims=(0, 250), ylims=(0, 250), xlabel="[cI] (nM)", ylabel="[cro] (nM)", aspect_ratio=:equal, legend=:top, size=(600, 600))

# ## Fig 7.11 (B)
ps711b = ComponentArray(ps711; delta_r = 0.2)
prob2 = remake(prob, p = ps711b)

rx = range(0, 250, 251)
ry = range(0, 250, 251)
rxy = range(0, 250, 21)
xx = [x for y in rxy, x in rxy]
yy = [y for y in rxy, x in rxy]

∂R2 = (x, y) -> model711!(ComponentArray(r=0.0, c=0.0), ComponentArray(r=x, c=y), ps711b, nothing).dr
∂C2 = (x, y) -> model711!(ComponentArray(r=0.0, c=0.0), ComponentArray(r=x, c=y), ps711b, nothing).dc
∂F2 = function (x, y; scale=0.2)
    dR, dC = model711!(ComponentArray(r=0.0, c=0.0), ComponentArray(r=x, c=y), ps711b, nothing)
    s = sqrt(hypot(dR, dC)) * scale
    return (dR / s, dC / s)
end

fig = plot(title="Figure 7.11 (B)")
contour!(fig, rx, ry, ∂R2, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="R nullcline")
contour!(fig, rx, ry, ∂C2, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="C nullcline")
quiver!(fig, xx, yy, quiver=∂F2, line=(:lightblue), arrow=(:closed))

plot!(fig, xlims=(0, 250), ylims=(0, 250), xlabel="[cI] (nM)", ylabel="[cro] (nM)", aspect_ratio=:equal, legend=:top, size=(600, 600))
