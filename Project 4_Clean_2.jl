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
    
    if var == theta
    varopt=optimizetheta(;thetaopt = varopt)
    return varopt
    end

    if var == x
        varopt=optimizex(;xopt = varopt)
        return varopt
    end

end

# if thetaopt !=0

# thetaopt =run(;varopt=thetaopt, var = theta)

# else

#     thetaopt =run(var = theta)

# end

if xopt !=0

    xopt =run(;varopt=xopt, var = x)

else

    xopt =run(var = x)

end