module Startup

using PrecompileTools: @recompile_invalidations, @setup_workload, @compile_workload
using ModelingToolkit
using Catalyst
using OrdinaryDiffEq
using Plots

@setup_workload begin
    function collins_sys(; name=:collins)
        @independent_variables t
        D = Differential(t)
        @parameters a1=3 a2=2.5 β=4 γ=4
        @variables i1(t) i2(t) s1(t)=0.075 s2(t)=2.5
        hil(x, k) = x / (x + k)
        hil(x, k, n) = hil(x^n, k^n)
        eqs = [
            D(s1) ~ a1 * hil(1 + i2, s2, β) - s1,
            D(s2) ~ a2 * hil(1 + i1, s1, γ) - s2,
            i1 ~ 10 * (t > 30) * (t < 40),
            i2 ~ 10 * (t > 10) * (t < 20)
        ]
        return System(eqs, t; name)
    end

    up303 = Dict(:S => 5.0, :ES => 0.0, :P => 0.0, :E=>1.0, :k1 => 30.0, :km1 => 1.0, :k2 => 10.0)

    @compile_workload begin
        @mtkbuild sys = collins_sys()
        prob = ODEProblem(sys, [], 50.0)
        sol = solve(prob, TRBDF2(), tstops=[10.0, 20.0, 30.0, 40.0])
        plot(sol, title="Fig. 1.7", xlabel="Time", ylabel="Concentration")

        rn303 = @reaction_network begin
            k1, S + E --> ES
            km1, ES --> S + E
            k2, ES --> P + E
        end
        prob303 = ODEProblem(rn303, up303, (0.0, 1.0); remove_conserved=true)
        sol303 = solve(prob303, TRBDF2())
    end
end

end # module Startup
