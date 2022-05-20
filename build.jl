# Adapted from https://github.com/terasakisatoshi/sysimage_creator/
import Pkg

Pkg.add(["PackageCompiler", "IJulia"])

using PackageCompiler

sysimage_path = joinpath(@__DIR__, "sysimage.so")

@info "SysImage path: " sysimage_path

# Do not include too many packages because memory on CI machines is limited
PackageCompiler.create_sysimage(
    ["Plots"];
    project=".",
    sysimage_path=sysimage_path,
    cpu_target=PackageCompiler.default_app_cpu_target()
)

nthreads = Sys.CPU_THREADS

using IJulia
kernelpath = IJulia.installkernel("Julia-sys-$(nthreads)-threads",
    "--project=@.", "--sysimage=$(sysimage_path)";
    specname="julia", env=Dict("JULIA_NUM_THREADS" => "$(nthreads)",)
)

@info "IJulia kernal path: " kernelpath

Pkg.rm("PackageCompiler")
Pkg.gc()
