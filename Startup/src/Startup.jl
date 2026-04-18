module Startup

using PrecompileTools: @recompile_invalidations

@recompile_invalidations begin
    using ModelingToolkit
    using OrdinaryDiffEq
    using SteadyStateDiffEq
    using Catalyst
    using Plots
    using CSV
    using DataFrames
end

end # module Startup
