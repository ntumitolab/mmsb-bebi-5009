# Adapted from https://github.com/terasakisatoshi/sysimage_creator/
using Pkg
Pkg.instantiate()

using IJulia
nthreads = Threads.nthreads()

IJulia.installkernel(
    "Julia",
    "--project=@.";
    env=Dict("JULIA_NUM_THREADS" => "$(nthreads)",)
)
