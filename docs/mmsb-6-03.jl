#===
# Fig 6.3

Two component pathway
===#

using ModelingToolkit
using Catalyst
using DifferentialEquations
using Plots
Plots.default(linewidth=2)

#---

function init_model63(events=[]; name)
    @variables t
    vs = @variables R(t) P(t) RL(t) Ps(t)
    ps = @parameters Rt=3 Pt=8 k1=5 km1=1 k2=6 k3=2 L=0
    D = Differential(t)
    eqs = [
        Rt ~ R + RL,
        Pt ~ P + Ps,
        D(R) ~ -k1 * L * R + km1 * RL,
        D(P) ~ -k2 * P * RL + k3 * Ps
    ]
    sys = ODESystem(eqs, t, vs, ps; name, discrete_events=events)
    return structural_simplify(sys)
end

# Fig. 6.3 A
@variables t
@parameters L
add_ligand = [1.0] => [L ~ 3.0]
rm_ligand = [3.0] => [L ~ 0.0]

@named sys = init_model63([add_ligand, rm_ligand])
@unpack R, Rt, P, Pt = sys
u0 = [R => Rt, P => Pt]
tend = 10.
prob = ODEProblem(sys, u0, tend)

#---

sol = solve(prob)

@unpack RL, Ps = sys
fig = plot(
    sol, idxs=[RL, Ps],
    labels= ["RL" "P*"],
    title="Fig. 6.3 (A)",
    xlabel="Time",
    ylabel="Concentration"
)

plot!(fig, t -> 3 * (1<=t<=3), label="Ligand", line=(:black, :dash), linealpha=0.7)

# Fig 6.3 (B)

@named sys = init_model63()
@unpack L, Ps, RL = sys
idx = findfirst(isequal(L), parameters(sys))
lrange = 0:0.01:1

prob_func = function(prob, i, repeat)
    p = copy(prob.p)
    p[idx] = lrange[i]
    remake(prob, p=p)
end

#---
prob = SteadyStateProblem(sys, u0, [])
#---
eprob = EnsembleProblem(prob; prob_func)
#---
sim = solve(eprob, DynamicSS(Rodas5()); trajectories=length(lrange));

#---

pstar = map(s->s[Ps], sim)
rl = map(s->s[RL], sim)
plot(lrange, [pstar rl], label=["P*" "RL"], title="Fig. 6.3 (B)",
xlabel="Ligand", ylabel="Steady-state concentration", xlims=(0, 1), ylims=(0, 8))
