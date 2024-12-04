using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
include("OptimizeChordsFunctions.jl")

function OptimizeChord(num_sec, scale_factor)

    chord_opt, fopt, info = RunOptimizer(num_sec, scale_factor)

    chords, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, scale_factor, chord_opt)

#    chords = zeros(num_sec+1)

#    for i in 1:num_sec+1
#    chords[i] = chord_opt[i]
#    end

#    Vinf=5*chord_opt[num_sec+2]

    xle_opt = XLEcalc(xle, chords, num_sec)

    symmetric = true

    Sref, ref = ReferenceCalculation(num_sec, yle, span, Vinf, chords)

    fs = FreestreamParams(Vinf)

    # reconstruct surface
    grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chords, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces_opt = [surface_opt]

    system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

    properties_opt = get_surface_properties(system_opt)

    write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=true)

    println("Optimized leading edge values: ", xle_opt)
    println("Optimized chord values: ", chords)
    println("Optimized Vinf = ", Vinf)

    # Plot the chords
    plot_chords(xle_opt, yle, chord_opt, num_sec)
    savefig("Chord Plot")
end

OptimizeChord(12, 55)