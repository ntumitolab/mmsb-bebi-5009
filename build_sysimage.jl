import Pkg
Pkg.add(["PackageCompiler", "IJulia"])
using PackageCompiler
PackageCompiler.create_sysimage(
    ["Plots", "DifferentialEquations", "ModelingToolkit", "Catalyst"];
    project=".", sysimage_path="$(pwd())/sysimage.so")

using IJulia
IJulia.installkernel("Julia", "--sysimage=$(pwd())/sysimage.so")

Pkg.rm("PackageCompiler")
Pkg.gc()
