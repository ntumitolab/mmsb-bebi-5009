# # Fig 6.14
# Model of E. coli chemotaxis signalling pathway
using OrdinaryDiffEq
using Catalyst
using Plots

@time "Build system" rn614 = @reaction_network begin
    @discrete_events begin
        (t == 10.0) => [L => 40.0]
        (t == 30.0) => [L => 80.0]
    end
    mm(Am, k1 * BP, KM1), Am => A
    mm(AmL, k2 * BP, KM2), AmL => AL
    km1 * R , A => Am
    km2 * R , AL => AmL
    (k3 * L, km3), Am <--> AmL
    (k4 * L, km4), A <--> AL
    (k5 * Am, km5), B <--> BP
end

#---
alg = FBDF()
up614 = Dict(
    :Am => 0.0360,
    :AmL => 1.5593,
    :A => 0.0595,
    :AL => 0.3504,
    :B => 0.7356,
    :BP => 0.2644,
    :k1 => 200.0,
    :k2 => 1.0,
    :k3 => 1.0,
    :km1 => 1.0,
    :km2 => 1.0,
    :km3 => 1.0,
    :k4 => 1.0,
    :km4 => 1.0,
    :k5 => 0.05,
    :km5 => 0.005,
    :KM1 => 1.0,
    :KM2 => 1.0,
    :L => 20.0,
    :R => 1.0,
)

tend = 50.0
@time "Build problem" prob614 = ODEProblem(rn614, up614, (0.0, tend); remove_conserved=true)
@time "Solve problem" sol614 = solve(prob614, alg; tstops=[10.0, 30.0])
plot(sol614, idxs=[:Am], title="Fig 6.14", xlabel="Time", ylabel="Concentration")
