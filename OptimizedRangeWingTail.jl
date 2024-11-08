using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
using LinearAlgebra
include("OptimizedRangeWingTailFunctions.jl")

function OptimizeChord(num_sec, scale_factor, xstart)

    xopt, fopt, info = RunOptimizer(num_sec, scale_factor, xstart)

    dt, chord_opt, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, scale_factor, xopt)

    xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v=SetUpTail()

    Sref, ref = ReferenceCalculation(num_sec, yle, span, Vinf, chord_opt)

    fs = FreestreamParams(Vinf)

    # construct surface
    grid_opt, surface_opt = wing_to_surface_panels(xle, yle, zle, chord_opt, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid_opt, htail_opt = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, [dt, 0.0, 0.0])
    translate!(htail, [dt, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid_opt, vtail_opt = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vgrid, [dt, 0.0, 0.0])
    translate!(vtail, [dt, 0.0, 0.0])

            # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = [true, true, false]

    # create vector containing all surfaces
    surfaces_opt = [surface_opt, htail_opt, vtail_opt]

    system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

    properties_opt = get_surface_properties(system_opt)

    write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

    println("Optimized leading edge values: ", xle)
    println("Optimized chord values: ", chord_opt)
    println("Optimized flight speed:  ", Vinf)
    println("Optimized tail distance:  ", dt)

    # Plot the chords
    plot_chords(xle, yle, chord_opt, num_sec)
    savefig("Chord Plot")
    return xopt
end

N=4
if @isdefined(xstart)
    if length(xstart) == N+3
    I=xstart
    xstart=OptimizeChord(N, 1, I)
    elseif length(xstart) > N+3
        I=zeros(N+3)
        for i in 1:N+3
        I[i]=xstart[i]
        end
        xstart=OptimizeChord(N, 1, I)
    else
        I=zeros(N+3)
        for i in 1:length(xstart)
        I[i]=xstart[i]
        end
        xstart=OptimizeChord(N, 1, I)
    end
else
    I=zeros(N+3)
    xstart=OptimizeChord(N, 1, I)
end