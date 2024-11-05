using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
include("Project5Functions.jl")

function OptimizeWing(num_sec, scale_factor)

    x_opt, fopt, info = RunOptimizer(num_sec, scale_factor)

    for i in 1:num_sec+1
        theta_opt[i] = x_opt[i+num_sec+1]/19
    end

    span_opt = 4 * xopt[num_sec+3]

    for i in 1:num_sec+1
        chord_opt[i] = x_opt[i]
    end

    span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, scale_factor)

    xle_opt = XLEcalc(xle, chord_opt, num_sec)

    symmetric = true

    Sref, ref = ReferenceCalculation(num_sec, yle, span_opt, Vinf, chord_opt, theta_opt)

    fs = FreestreamParams(Vinf)

    # reconstruct surface
    grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chord_opt, theta_opt, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces_opt = [surface_opt]

    system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

    properties_opt = get_surface_properties(system_opt)

    write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=true)

    println("Optimized leading edge values: ", xle_opt)
    println("Optimized chord values: ", chord_opt)

    # Plot the chords
    plot_chords(xle_opt, yle, chord_opt, num_sec)
    savefig("Chord Plot")
end

OptimizeChord(12, 5)