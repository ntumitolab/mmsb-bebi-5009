md"""
# Plotting with Julia

- [Plots.jl](https://github.com/JuliaPlots/Plots.jl): powerful and convenient visualization mwith ultiple backends. See also [Plots.jl docs](https://docs.juliaplots.org/latest/)
- [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl): `matplotlib` in Julia. See also [matplotlib docs](https://matplotlib.org/stable/index.html)
- [Makie.jl](https://github.com/MakieOrg/Makie.jl): a data visualization ecosystem for the Julia programming language, with high performance and extensibility. See also [Makie.jl docs](https://docs.makie.org/stable/)
"""

using Plots

f(x) = sin(sin(x) + 1)

# Prepare data
xs = 0.0:0.1:4pi
ys = f.(xs)

# Line plots connect the data points
plot(xs, ys)

# scatter plots show the data points only
scatter(xs, ys)

# you can trace functions directly
plot(f, xs)

# Trace a function within a range
plot(f, 0.0, 4pi)

# Customization example
plot(f, xs,
     label="My line", legend=:bottom,
     title="My Title",  line=(:red, 3),
     xlim = (0.0, 5.0), ylim = (-1.0, 1.5),
     xlabel="Time", ylabel="My Mood", border=:box)

# Multiple series: each row is one observation; each column is a variable.
f2(x) = cos(cos(x) + 1)
y2 = f2.(xs)
plot(xs, [ys y2])

# Plotting two functions with customizations
plot(xs, [f, f2], label=["f1" "f2"], linecolor=[:black :green], title="Two time series")

# Building the plot in multiple steps in the object-oriented way
xMin = 0.0
xMax = 4.0π
p1 = plot(f, xMin, xMax, label="f1", lc=:black)
plot!(p1 , f2, xMin, xMax, label="f2", lc=:lightsalmon)
plot!(p1, title = "My title", legend=:outertop)

# Parametric plot
xₜ(t) = sin(t)
yₜ(t) = sin(2t)

plot(xₜ, yₜ, 0, 2π, leg=false, fill=(0,:orange))

# Subplots
p1 = plot(f, xs)
p2 = plot(f2, xs)
plot(p1, p2)

# Subplot layout
psub = plot(p1, p2, layout=(2, 1))

md"""
## Vector field

### `Plots.jl`

```julia
# Quiver plot
quiver(vec(x2d), vec(y2d), quiver=(vec(vx2d), vec(vy2d))
# Or if you have a function f(x,y) -> (vx, vy)
quiver(x2d, y2d, quiver=f)
```

### `PyPlot.jl`:

```julia
using PyPlot as plt
plt.quiver(X2d, Y2d, U2d, V2d)
```

See also: [matplotlib: quiver plot](https://matplotlib.org/stable/gallery/images_contours_and_fields/quiver_demo.html#sphx-glr-gallery-images-contours-and-fields-quiver-demo-py)
"""

using Plots

## ∇ = \nabla <TAB>
function ∇f(x, y; scale=hypot(x, y)^0.5 * 3)
    return [-y, x] ./ scale
end

## x and y grid points
r = -1.0:0.2:1.0
xx = [x for y in r, x in r]
yy = [y for y in r, x in r]

quiver(xx, yy, quiver=∇f, aspect_ratio=:equal, line=(:black), arrow=(:closed))

#===
## Save figure

`savefig(filename)`
===#

savefig("vfield.png")
savefig(psub, "subplots.png")

# ## Runtime information

import InteractiveUtils
InteractiveUtils.versioninfo()

#---

import Pkg
Pkg.status()
