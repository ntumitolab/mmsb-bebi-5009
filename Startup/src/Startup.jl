module Startup

import PrecompileTools
PrecompileTools.@recompile_invalidations begin
    using ModelingToolkit
    using Catalyst
    using OrdinaryDiffEq
    using Plots
    using DataFrames
    using CSV
    using SteadyStateDiffEq
end

end # module Startup
