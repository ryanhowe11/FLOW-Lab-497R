using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
using Statistics
using LinearAlgebra
include("OptimizedRangeWingTwistTailDihedralSpanFunctions.jl")

function OptimizeTwist(num_sec, density, xstart)

    xopt, fopt, info = RunOptimizer(num_sec, density, xstart)

    lh, lv, dt, chord_opt, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, density, xopt)

    xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v=SetUpTail(xopt, lh, lv)

    Sref, ref = ReferenceCalculation(num_sec, yle, span, Vinf, chord_opt)

    println(Sref)

    fs = FreestreamParams(Vinf)

    # construct surface
    grid_opt, surface_opt = wing_to_surface_panels(xle, yle, zle, chord_opt, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid_opt, htail_opt = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    VortexLattice.translate!(hgrid_opt, [dt, 0.0, 0.0])
    VortexLattice.translate!(htail_opt, [dt, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid_opt, vtail_opt = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    VortexLattice.translate!(vgrid_opt, [dt, 0.0, 0.0])
    VortexLattice.translate!(vtail_opt, [dt, 0.0, 0.0])

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
    println("Optimized Horizontal tail size = ", lh)
    println("Optimized Vertical tail size = ", lv)
    # Plot the chords
    plot, plot_chords(xle, yle, chord_opt, num_sec)
    savefig("Chord Plot")
    return xopt, fopt, plot
end

function GetVTKandValues(num_sec, density, xstart)

    lh, lv, dt, chord_opt, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, density, xstart)

    xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v=SetUpTail(xstart, lh, lv)

    Sref, ref = ReferenceCalculation(num_sec, yle, span, Vinf, chord_opt)

    println(Sref)

    fs = FreestreamParams(Vinf)

    # construct surface
    grid_opt, surface_opt = wing_to_surface_panels(xle, yle, zle, chord_opt, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid_opt, htail_opt = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    VortexLattice.translate!(hgrid_opt, [dt, 0.0, 0.0])
    VortexLattice.translate!(htail_opt, [dt, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid_opt, vtail_opt = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    VortexLattice.translate!(vgrid_opt, [dt, 0.0, 0.0])
    VortexLattice.translate!(vtail_opt, [dt, 0.0, 0.0])

            # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = [true, true, false]

    # create vector containing all surfaces
    surfaces_opt = [surface_opt, htail_opt, vtail_opt]

    system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

    properties_opt = get_surface_properties(system_opt)

    write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

    println("Optimized leading edge values: ", xle)
    println("Optimized chord values: ", chord_opt)
    println("Optimized twist value:  ", theta*180/pi)
    println("Optimized dihedral angle:  ", phi*180/pi)
    println("Optimized flight speed:  ", Vinf)
    println("Optimized tail distance:  ", dt)
    println("Optimized Horizontal tail size = ", lh)
    println("Optimized Vertical tail size = ", lv)
    # Plot the chords
    plot=plot_chords(xle, yle, chord_opt, num_sec)
    savefig("Chord Plot")
    return surfaces_opt, properties_opt, symmetric
end

density=50
N=10

# if @isdefined(xstart)
#     if length(xstart) == N+5
#     xstart_vec=xstart
#     xstart=OptimizeChord(N, density, xstart_vec)
#     elseif length(xstart) > N+5
#         xstart_vec=ones(N+5)
#         for i in 1:N+5
#         xstart_vec[i]=xstart[i]
#         end
#         xstart=OptimizeChord(N, density, xstart_vec)
#     else
#         xstart_vec=ones(N+5)
#         for i in 1:length(xstart)
#             xstart_vec[i]=xstart[i]
#         end
#         xstart=OptimizeChord(N, density, xstart_vec)
#     end
#     else
#     xstart_vec=ones(N+5)
#     xstart=OptimizeChord(N, density, xstart_vec)
# end

xstart_vec = 0.2*ones(N+7)
x1, f1, p1 = OptimizeTwist(N, density, xstart_vec)
savefig(p1, "OWS_chords_plot_1.png")

xstart_vec = ones(N+7)
x2, f2, p2 = OptimizeTwist(N, density, xstart_vec)
savefig(p2, "OWS_chords_plot_2.png")

xstart_vec = 5*ones(N+7)
x3, f3, p3 = OptimizeTwist(N, density, xstart_vec)
savefig(p3, "OWS_chords_plot_3.png")

xstart_vec = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 2, 2, 2, 2, 2, 2, 5]
x4, f4, p4 = OptimizeTwist(N, density, xstart_vec)
savefig(p4, "OWS_chords_plot_4.png")

xstart_vec = [2.1, 1.9, 1.7, 1.5, 1.3, 1.1, .9, .7, .5, .3, 2, 2, 2, 2, 2, 2, 5]
x5, f5, p5 = OptimizeTwist(N, density, xstart_vec)
savefig(p5, "OWS_chords_plot_5.png")

xstart_vec = [2.1, 2.9, 1.7, 3.5, .3, 4.1, 2.9, 1.7, 3.5, 4.73, 2, 2, 2, 2, 2, 2, 5]
x6, f6, p6 = OptimizeTwist(N, density, xstart_vec)
savefig(p6, "OWS_chords_plot_6.png")

xstart_vec = [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1, 2, 3, 4, 5, 6, 7]
x7, f7, p7 = OptimizeTwist(N, density, xstart_vec)
savefig(p7, "OWS_chords_plot_7.png")

xstart_vec = [2.1, 1.9, 1.7, 1.5, 1.3, 1.1, .9, .7, .5, .3, 6, 5, 4, 3, 2, 1, 3]
x8, f8, p8 = OptimizeTwist(N, density, xstart_vec)
savefig(p8, "OWS_chords_plot_8.png")

xstart_vec = [2.1, 2.9, 1.7, 3.5, .3, 4.1, 2.9, 1.7, 3.5, 4.73, 3, 2, 4, 1, 5, 2, 3]
x9, f9, p9 = OptimizeTwist(N, density, xstart_vec)
savefig(p9, "OWS_chords_plot_9.png")

xstart_vec = [1.1, 3.9, .07, 4.56, .03, 4.91, 3.9, 2.1, .5, 4.73, .1, 4.8, 2.5, 6, 1, 5, 7]
x10, f10, p10 = OptimizeTwist(N, density, xstart_vec)
savefig(p10, "OWS_chords_plot_10.png")

# xstart_vec = x1
# x1, f1, p1 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p1, "OWS_chords_plot_1.png")

# xstart_vec = x2
# x2, f2, p2 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p2, "OWS_chords_plot_2.png")

# xstart_vec = x4
# x4, f4, p4 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p4, "OWS_chords_plot_4.png")

# xstart_vec = x5
# x5, f5, p5 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p5, "OWS_chords_plot_5.png")

# xstart_vec = x6
# x6, f6, p6 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p6, "OWS_chords_plot_6.png")

# surfaces_opt, properties_opt, symmetric = GetVTKandValues(N, density, x4)
# write_vtk("Plane_4", surfaces_opt, properties_opt; symmetric=symmetric)

# surfaces_opt, properties_opt, symmetric = GetVTKandValues(N, density, x5)
# write_vtk("Plane_5", surfaces_opt, properties_opt; symmetric=symmetric)

# surfaces_opt, properties_opt, symmetric = GetVTKandValues(N, density, x6)
# write_vtk("Plane_6", surfaces_opt, properties_opt; symmetric=symmetric)

# xstart_vec = x7
# x7, f7, p7 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p7, "OWS_chords_plot_7.png")

# xstart_vec = x8
# x8, f8, p8 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p8, "OWS_chords_plot_8.png")

# xstart_vec = x9
# x9, f9, p9 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p9, "OWS_chords_plot_9.png")

# xstart_vec = x10
# x10, f10, p10 = OptimizeChordTwist(N, density, xstart_vec)
# savefig(p10, "OWS_chords_plot_10.png")