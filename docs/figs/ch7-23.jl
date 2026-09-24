# # Fig 7.23 and 7.31
# Hasty synthetic oscillator model
using Startup: meshgrid, get_gradient, normalize_gradient
using Catalyst
using DifferentialEquations
using Plots
Plots.gr(linewidth=1.5)

#---
v0_723(x, y, alpha, sigma) = (1 + x^2 + alpha * sigma * x^4) / ((1 + x^2 + sigma * x^4) * (1 + y^4))
@time "Build system" rn723 = @reaction_network begin
    v0_723(x, y, alpha, sigma), 0 --> x
    ay * v0_723(x, y, alpha, sigma), 0 --> y
    gammax, x --> 0
    gammay, y --> 0
end

# ## Fig. 7.23 (A)
alg = FBDF()
up723 = Dict(
    :alpha => 11.0, :sigma => 2.0, :gammax => 0.2, :gammay => 0.012, :ay => 0.2,
    :x => 0.3963, :y => 2.3346,
)

tend = 300.0
@time "Build problem" prob723 = ODEProblem(rn723, up723, (0.0, tend))
@time "Solve problem" sol723 = solve(prob723, FBDF())
plot(sol723, idxs=[:x, :y], xlabel="Time", ylabel="Concentration", title="Fig 7.23 (A)", legend=:right)

# ## Fig. 7.23 (B)
# Vector field
xrange = range(0, 1.5, 51)
yrange = range(0, 3, 51)
@unpack x, y = rn723

xx, yy = meshgrid(xrange, yrange)
dx, dy = get_gradient(prob723, x, y, xx, yy)
contour(xrange, yrange, dx, levels=[0], line=(:black, :solid), colorbar=false)
plot!([], [], line=(:black, :solid), label="x nullcline")
contour!(xrange, yrange, dy, levels=[0], line=(:black, :dash), colorbar=false)
plot!([], [], line=(:black, :dash), label="y nullcline")
plot!(sol723, idxs=(x, y), label="Trajectory")
plot!(xlabel="x", ylabel="y", title="Fig 7.23 (B)", legend=:bottomright, size=(600, 600))

# ## Fig. 7.31
@time "Build system" rn731 = @reaction_network begin
    v0_723(x1, y1, alpha, sigma), 0 --> x1
    ay * v0_723(x1, y1, alpha, sigma), 0 --> y1
    v0_723(x2, y2, alpha, sigma), 0 --> x2
    ay * v0_723(x2, y2, alpha, sigma), 0 --> y2
    gammax, x1 --> 0
    gammay, y1 --> 0
    gammax, x2 --> 0
    gammay, y2 --> 0
    (D, D), x1 <--> x2
    (D, D), y1 <--> y2
end

#---
alg = FBDF()
up731 = Dict(
    :alpha => 11.0,
    :sigma => 2.0,
    :gammax => 0.2,
    :gammay => 0.012,
    :ay => 0.2,
    :D => 0.015,
    :x1 => 0.3963,
    :y1 => 2.3346,
    :x2 => 0.5578,
    :y2 => 1.9317,
)

tend = 500.0
@time "Build problem" prob731 = ODEProblem(rn731, up731, (0.0, tend))
@time "Solve problem" sol731 = solve(prob731, alg)
plot(sol731, xlabel="Time (a.u.)", ylabel="Concentration (a.u.)", title="Fig. 7.31", legend=:left)
