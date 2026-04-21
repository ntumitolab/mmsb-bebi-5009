module Startup

import PrecompileTools
PrecompileTools.@recompile_invalidations begin
    using Catalyst
    using CSV
    using DataFrames
    using DataInterpolations
    using DiffEqCallbacks
    using ForwardDiff
    using JumpProcesses
    using Latexify
    using ModelingToolkit
    using OrdinaryDiffEq
    using Plots
    using SteadyStateDiffEq
end

end # module Startup
