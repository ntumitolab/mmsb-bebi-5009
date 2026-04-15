
using OrdinaryDiffEq
using Plots
using Catalyst
using Latexify

glycolysis = @reaction_network begin
    hillr(ATP, k1, Ki1A, 4) * mm(ATP, 1, K1A), ATP => G6P + ADP
    k2, G6P --> 0
    k3 * mm(F6P, 1, K3F) * mm(ATP, 1, K3A), F6P + ATP => F26P2 + ADP
    (k4, k4), G6P <--> F6P
    k5 * mm(F6P, 1, K5F) * mm(ATP, 1, K5A) * (1 + mm(F26P2, alphaB, K5B)) * hill(ADP, 1, K5D, 4) * hillr(ATP, 1, Ki5A, 4), F6P + ATP => 2TP + ADP
    k6 * mm(TP, 1, K6T) * hill(ADP, 1, K6D, 2), TP + 2ADP => 2ATP
    (k7f, k7r), 2ADP <--> ATP + AMP
    k8, F26P2 --> F6P
    k9, ATP --> ADP
end

unknowns(glycolysis)
netstoichmat(glycolysis)
conservedequations(glycolysis)
println(latexify(glycolysis))

up = [:k1=>1.0, :Ki1A=>2.0, :K1A=>0.1, :k2=>0.05, :k3=>0.005, :K3F=>0.15, :K3A=>0.25, :k4=>1.0, :k5=>1.5, :K5F=>0.05, :K5A=>0.05, :alphaB=>4.0, :K5B=>0.03, :K5D=>0.4, :Ki5A=>1.0, :k6=>1.0, :K6T=>0.5, :K6D=>0.2, :k7f=>10.0, :k7r=>10.0, :k8=>0.02, :k9=>0.2, :G6P=>0.0, :F6P=>0.0, :F26P2=>0.0, :TP=>0.0, :ATP=>1.8, :ADP=>0.6, :AMP=>0.2]

tend = 1000.0
prob = ODEProblem(glycolysis, up, (0.0, tend); combinatoric_ratelaw=false)
@time sol = solve(prob, TRBDF2(); abstol=1e-8, reltol=1e-6)

p1 = plot(sol, idxs=[:ATP, :ADP, :AMP], xlabel="Time", ylabel="Concentration", linewidth=2, title="Energy states", )

p2 = plot(sol, idxs=[:G6P, :F6P, :F26P2, :TP], xlabel="Time", ylabel="Concentration", linewidth=2, title="Energy states")

plot(p1, p2, layout=(2, 1), legend=:topright, size=(600, 600))
