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

function run()
    num_sec = 24
    sec_2 = Int(.5*num_sec)
    scale_factor = 1
    Prev_Run = 1
    Include_Airframe = 1

    #Creating the optimization problem
    function wing_optimizer(g, theta)

        thetaopt=theta

        D, weight, rho, Vinf, Sref, CL, chords = Drag_Calculation_Theta(thetaopt)

        g[1]=weight-.5*rho*Vinf^2*Sref*CL

        # # Calculate chord differences
        # for i in 1:num_sec
        #     g[i+1] = thetaopt[i] - thetaopt[i+1]#+.25/num_sec
        # end

    return D
    end

    theta0, ng, ltheta, utheta, lg, ug, options, g = setup(Prev_Run)

    thetaopt, fopt, info = minimize(wing_optimizer, theta0, ng, ltheta, utheta, lg, ug, options)

    println("Optimized twist values: ", thetaopt*180/pi)

    if Include_Airframe == 1
        
        span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, dh, dv = AirFrameVariables(scale_factor)
    
    else

        span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Variables(scale_factor)
    
    end

    theta=thetaopt

    #Reference Parameters
    Sref, cref, rref = GetRef(c, chords, theta, span)

    ref = Reference(Sref, cref, span, rref, Vinf)

    xle_opt = XLE_calc(chords, num_sec)

    if Include_Airframe == 1

        # reconstruct surface
        wgrid_opt, wing_opt = wing_to_surface_panels(xle_opt, yle, zle, chords, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

        # generate surface panels for horizontal tail
        hgrid_opt, htail_opt = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
        VortexLattice.translate!(hgrid_opt, [dh, 0.0, 0.0])
        VortexLattice.translate!(htail_opt, [dh, 0.0, 0.0])

        # generate surface panels for vertical tail
        vgrid_opt, vtail_opt = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
        VortexLattice.translate!(vgrid_opt, [dv, 0.0, 0.0])
        VortexLattice.translate!(vtail_opt, [dv, 0.0, 0.0])

        # now set normal vectors manually
        ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

        # overwrite normal vector for each wing panel
        for i = 1:length(wing_opt)
        wing_opt[i] = set_normal(wing_opt[i], ncp)
        end

        grids = [wgrid_opt, hgrid_opt, vgrid_opt]
        # create vector containing all surfaces
        surfaces_opt = [wing_opt, htail_opt, vtail_opt]
        surface_idopt = [1, 2, 3]
        
        system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric, surface_id=surface_idopt)
        
        properties_opt = get_surface_properties(system_opt)
        
        r, c = lifting_line_geometry(grids, 0.25)
        
        cf, cm = lifting_line_coefficients(system_opt, r, c; frame=Wind())
        
        symmetric = [true, true, false]
        
        write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

        CF, CM = body_forces(system_opt; frame=Stability())
        CDiff = far_field_drag(system_opt)
        CD, CY, CL = CF
        Cl, Cm, Cn = CM

        # Assuming cf is your vector of matrices
        z_direction_coefficients = []

        # Loop through each matrix in the cf vector
        for matrix in cf
            # Extract the third row (z-direction force coefficients)
            z_coefficients = matrix[3, :]
            push!(z_direction_coefficients, z_coefficients)
        end

        y2 = collect(range(0, stop=span, step=0.1))

        # Convert the list of z-direction coefficients to an array if needed
        lifting_coefficients = hcat(z_direction_coefficients...)

        Lift_prime = .5*rho*Vinf^2*Chord_Length*lifting_coefficients[:, 1]
        push!(Lift_prime, 0)

        yle_2 = [i * (span / (num_sec)) for i in 0:(num_sec)]

        area_prime = trapz(yle_2, Lift_prime)

        bprime = span
        aprime = 4*area_prime/(bprime*pi)

        # Calculate the ideal elliptic lift distribution
        θ = range(0, π/2, length=100)
        x = bprime * cos.(θ)
        y = aprime * sin.(θ)
        Cl_max = maximum(Lift_prime)
        elliptical_distribution = Cl_max * sqrt.(1 .- (y2 ./ span).^2)

        # Plot the lift distribution
        plot(yle_2, Lift_prime, label="Optimized Lift Distribution", xlabel="Spanwise Location (y)", ylabel="Lift Coefficient (Cl)")
        plot!(y2, elliptical_distribution, label="Elliptical Lift Distribution", linestyle=:dash)

        # Save the lift distribution plot as a PDF
        savefig("Lift_Distribution_along_the_Span_Twist_Optimization_With_Tail.pdf")

    else
        # reconstruct surface
        grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chords, thetaopt, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

        surfaces_opt, system_opt, properties_opt = VTK_setup(surface_opt, Vinf, alpha_angle, beta, Omega, ref, symmetric)

        write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

        frame=Wind()

        cf = LocalizedCoefficients(system_opt, grid_opt, frame)
        
        yle_2, y2, Lift_prime, elliptical_distribution = GetLiftCurves(cf, chords, rho, Vinf, num_sec, span)
    
        PlotLiftDistr(yle_2, y2, Lift_prime, elliptical_distribution)    

    end

end

run()