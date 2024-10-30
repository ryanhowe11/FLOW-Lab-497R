using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
using Calculus
using Trapz
using LinearAlgebra
include("Project4_functions.jl")
include("Project4_OptimizationFunctions.jl")

function run(;varopt, var)
    
    if var == "theta"
    varopt=optimizetheta(;thetaopt = varopt)
    return varopt
    end

    if var == "chord"
        varopt=optimizechord(;chordopt = varopt)
        return varopt
    end

end

# if thetaopt !=0

# thetaopt =run(;varopt=thetaopt, var = "theta")

# else

#     thetaopt =run(var = "theta")

# end


if chordopt !=0

    chordopt =run(;varopt=chordopt, var = "chord")

else

    chordopt =run(var = "chord")

end