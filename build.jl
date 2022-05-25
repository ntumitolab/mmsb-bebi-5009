# Adapted from https://github.com/terasakisatoshi/sysimage_creator/
import Pkg

Pkg.add(["IJulia"])
nthreads = Threads.nthreads()

using IJulia
kernelpath = IJulia.installkernel(
    "Julia",
    "--project=@.";
    env=Dict("JULIA_NUM_THREADS" => "$(nthreads)",)
)
