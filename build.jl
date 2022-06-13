# Adapted from https://github.com/terasakisatoshi/sysimage_creator/

nthreads = Threads.nthreads()

using IJulia
IJulia.installkernel(
    "Julia",
    "--project=@.";
    env=Dict("JULIA_NUM_THREADS" => "$(nthreads)",)
)
