module Startup

## PrecompileTools.jl tutorial: https://julialang.github.io/PrecompileTools.jl/stable/
using Reexport
@reexport using ModelingToolkit
@reexport using OrdinaryDiffEq
@reexport using Catalyst
using PrecompileTools: @setup_workload, @compile_workload
using ModelingToolkit: t_nounits as t, D_nounits as D

@compile_workload begin
    # Prepare data for `@compile_workload`
    # Collins toggle switch model
    function collins_sys(; name=:collins)
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
        return ODESystem(eqs, t; name)
    end

    # inside here, put a "toy example" of everything you want to be fast
    @mtkbuild sys = collins_sys()
    prob = ODEProblem(sys, [], (0.0, 50.0))
    sol = solve(prob, TRBDF2(), tstops=[10.0, 20.0, 30.0, 40.0])

    # Enzyme kinetics full model
    up303 = Dict(:S => 5.0, :ES => 0.0, :P => 0.0, :E=>1.0, :k1 => 30.0, :km1 => 1.0, :k2 => 10.0)

    rn303 = @reaction_network begin
        k1, S + E --> ES
        km1, ES --> S + E
        k2, ES --> P + E
    end
    prob303 = ODEProblem(rn303, up303, (0.0, 1.0); remove_conserved=true)
    sol303 = solve(prob303, TRBDF2())
end


end # module Startup
