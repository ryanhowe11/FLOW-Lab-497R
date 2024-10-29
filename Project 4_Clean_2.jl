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

function run(;thetaopt)
    num_sec = 12
    sec_2 = Int(.5*num_sec)
    Chord_Length = 1
    scale_factor = 1

    if thetaopt != 0
    Prev_Run = 1
    else
    Prev_Run = 0
    end
    
    Include_Airframe = 1
    ng = 13

    #Creating the optimization problem
    function wing_optimizer(g, theta)

        thetaopt=theta

        if Include_Airframe == 1
        
            D, weight, rho, Vinf, Sref, CL, chords, fs = Airframe_Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

        else

            D, weight, rho, Vinf, Sref, CL, chords, fs = Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

        end

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

    if Include_Airframe == 1
        
        span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, dh, dv = AirFrameVariables(scale_factor, num_sec, Chord_Length)
        D, weight, rho, Vinf, Sref, CL, chords, fs = Airframe_Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

    else

        span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Variables(scale_factor, num_sec, Chord_Length)
        D, weight, rho, Vinf, Sref, CL, chords, fs = Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

    end

    theta=thetaopt

    #Reference Parameters
    Sref, cref, rref = GetRef(c, chords, theta, span, num_sec, yle)

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

        grids = [wgrid_opt, hgrid_opt, vgrid_opt]
        # create vector containing all surfaces
        surfaces_opt = [wing_opt, htail_opt, vtail_opt]
        surface_idopt = [1, 2, 3]

        symmetric = [true, true, false]
        
        system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric, surface_id=surface_idopt)
        
        properties_opt = get_surface_properties(system_opt)
        
        r, c = lifting_line_geometry(grids, 0.25)
        
        cf, cm = lifting_line_coefficients(system_opt, r, c; frame=Wind())
        
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

        z_direction_coefficients_wing = z_direction_coefficients[1]

    y2 = collect(range(0, stop=span, step=0.1))
    # Convert the list of z-direction coefficients to an array if needed
    lifting_coefficients = hcat(z_direction_coefficients_wing...)
    lifting_coefficients = vec(transpose(lifting_coefficients))

    Lift_prime = .5*rho*Vinf^2 .*chords.*lifting_coefficients

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

    return thetaopt
end

if thetaopt !=0

thetaopt=run(;thetaopt=thetaopt)

else

    thetaopt=run()

end