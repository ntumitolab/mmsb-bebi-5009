import Pkg
Pkg.add(["PackageCompiler", "IJulia"])
using PackageCompiler
PackageCompiler.create_sysimage(
    ["Plots", "DifferentialEquations", "ModelingToolkit", "Catalyst"];
    project=".", sysimage_path="sysimage.so")

using IJulia
IJulia.installkernel("Julia", "--sysimage=sysimage.so")

Pkg.rm("PackageCompiler")
Pkg.gc()
