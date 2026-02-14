# # Fig 4.11
# Surface plots
using Plots
Plots.gr(framestyle = :box, linewidth=1.5)

#---
z1(x, y) = x^2 + 0.5y^2
z2(x, y) = (.2x^2-1)^2 + y^2
x1 = y1 = range(-1.0, 1.0, 41)
x2 = range(-2.75, 2.75, 41)
y2 = range(-0.75, 0.75, 41)
zz1 = [z1(x, y) for x in x1, y in y1]
zz2 = [z2(x, y) for x in x2, y in y2];

#---
pl1 = surface(x1, y1, zz1, title="Single-well potential", cbar=false)
pl2 = contourf(x1, y1, zz1, title="Single-well potential (contour)", cbar=false)
pl3 = surface(x2, y2, zz2, title="Double-well potential", cbar=false)
pl4 = contourf(x2, y2, zz2, title="Double-well potential (contour)", cbar=false)
plot(pl1, pl2, pl3, pl4, layout=(2, 2), size=(800, 800))
