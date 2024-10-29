using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
using Calculus
using Trapz
include("Project4_functions.jl")
function run()
    num_sec = 12
    sec_2 = Int(.5*num_sec)
    scale_factor = 1
    Chord_Length = 1
    Prev_Run=1
    ng = num_sec+1

    #Creating the optimization problem
    function wing_optimizer(g, theta)
        thetaopt=theta
        D, weight, rho, Vinf, Sref, CL, chords = Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)
        g[1]=weight-.5*rho*Vinf^2*Sref*CL
        # Calculate chord differences
        for i in 1:num_sec
            g[i+1] = abs(thetaopt[i]*180/pi) - abs(thetaopt[i+1]*180/pi)#+.25/num_sec
        end
    return D
    end
    theta0, ltheta, utheta, lg, ug, options, g = setup(Prev_Run, num_sec, ng)
    thetaopt, fopt, info = minimize(wing_optimizer, theta0, ng, ltheta, utheta, lg, ug, options)
    println("Optimized twist values: ", thetaopt*180/pi)
    span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Variables(scale_factor, num_sec, Chord_Length)
    theta=thetaopt
    #Reference Parameters
    Sref, cref, rref = GetRef(c, chords, theta, span, num_sec, yle)
    ref = Reference(Sref, cref, span, rref, Vinf)
    xle_opt = XLE_calc(chords, num_sec)
    # reconstruct surface
    grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chords, thetaopt, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)
    surfaces_opt, system_opt, properties_opt = VTK_setup(surface_opt, Vinf, alpha_angle, beta, Omega, ref, symmetric)
    write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)
    frame=Wind()
    cf = LocalizedCoefficients(system_opt, grid_opt, frame)
    yle_2, y2, Lift_prime, elliptical_distribution = GetLiftCurves(cf, chords, rho, Vinf, num_sec, span)
    PlotLiftDistr(yle_2, y2, Lift_prime, elliptical_distribution)
    return thetaopt
end
thetaopt=run()