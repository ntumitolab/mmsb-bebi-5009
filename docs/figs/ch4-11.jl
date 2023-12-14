#===
# Fig 4.11

Surface plots reference: [surface plots @ PlotsGallery.jl](https://goropikari.github.io/PlotsGallery.jl/src/surface.html)
===#

using Plots
Plots.default(linewidth=2)

# PNG output in Literate.jl
PNG(fig) = display("image/png", fig)

#---
z1(x, y) = x^2 + 0.5y^2
z2(x, y) = (.2x^2-1)^2 + y^2
x1 = y1 = range(-1.0, 1.0, 41)
x2 = range(-2.75, 2.75, 41)
y2 = range(-0.75, 0.75, 41)
p1 = surface(x1, y1, z1, title="Single-well potential")
p2 = contourf(x1, y1, z1)
p3 = surface(x2, y2, z2, title="Double-well potential")
p4 = contourf(x2, y2, z2)

fig = plot(p1, p2, p3, p4, size=(800, 600))

fig |> PNG
