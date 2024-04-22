#===
# Fig 7.11

Model of phage lambda decision switch
===#

using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---
function model711(u, p, t)
    r, c = u
    rd = r / 2 ## r Dimer
    cd = c / 2 ## c Dimer
    K1, K2, K3, K4 , delta_r, delta_c, a, b = p
    f1 = K1 * rd^2
    f2 = K2 * rd
    f3 = K3 * cd
    f4 = K4 * cd
    den = 1 + f1 * (1 + f2) + f3 * (1 + f4)
    dr = a * (1 + 10 * f1) / den - delta_r * r
    dc = b * (1 + f3) / den - delta_c * c
    return (dr, dc)
end

function model711!(D, u, p, t)
    D[1], D[2] = model711(u, p, t)
    return nothing
end

# ## Fig 7.11 (A)
ps1 = (K1=1, K2=0.1, K3=5, K4=0.5 , delta_r=0.02, delta_c=0.02, a=5, b=50)
tend = 6000.

rx = range(0, 250, 201)
ry = range(0, 250, 201)
rxy = range(0, 250, 21)
xx = [x for y in rxy, x in rxy]
yy = [y for y in rxy, x in rxy]

∂R = (x, y) -> model711((x, y), ps1, 0)[1]
∂C = (x, y) -> model711((x, y), ps1, 0)[2]

∂F = function (x, y; scale=0.2)
    dR, dC = model711((x, y), ps1, 0.0)
    s = sqrt(hypot(dR, dC)) * scale
    return (dR / s, dC / s)
end

fig = plot(title="Figure 7.11 (A)")
contour!(fig, rx, ry, ∂R, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="R nullcline")
contour!(fig, rx, ry, ∂C, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="C nullcline")
quiver!(fig, xx, yy, quiver=∂F, line=(:lightblue), arrow=(:closed))

plot!(fig, xlims=(0, 250), ylims=(0, 250), xlabel="[cI] (nM)", ylabel="[cro] (nM)", aspect_ratio=:equal, legend=:top, size=(600, 600))

# ## Fig 7.11 (B)
ps2 = merge(ps1, (; delta_r=ps1.delta_r * 10))

rx = range(0, 250, 201)
ry = range(0, 250, 201)
rxy = range(0, 250, 21)
xx = [x for y in rxy, x in rxy]
yy = [y for y in rxy, x in rxy]

∂R = (x, y) -> model711((x, y), ps2, 0)[1]
∂C = (x, y) -> model711((x, y), ps2, 0)[2]

∂F = function (x, y; scale=0.3)
    dR, dC = model711((x, y), ps2, 0.0)
    s = sqrt(hypot(dR, dC)) * scale
    return (dR / s, dC / s)
end

fig = plot(title="Figure 7.11 (B)")
contour!(fig, rx, ry, ∂R, levels=[0], cbar=false, line=(:black))
plot!(fig, Float64[], Float64[], line=(:black), label="R nullcline")
contour!(fig, rx, ry, ∂C, levels=[0], cbar=false, line=(:black, :dash))
plot!(fig, Float64[], Float64[], line=(:black, :dash), label="C nullcline")
quiver!(fig, xx, yy, quiver=∂F, line=(:lightblue), arrow=(:closed))

plot!(fig, xlims=(0, 250), ylims=(0, 250), xlabel="[cI] (nM)", ylabel="[cro] (nM)", aspect_ratio=:equal, legend=:top, size=(600, 600))
