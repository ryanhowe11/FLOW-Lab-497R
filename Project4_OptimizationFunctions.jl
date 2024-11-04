function optimizetheta(;thetaopt)

        num_sec = 12
        Chord_Length = 1
        scale_factor = 1

        if thetaopt != 0
        Prev_Run = 1
        else
        Prev_Run = 0
        end
        
        Include_Airframe = 0
        ng = num_sec+1

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

        theta0, ltheta, utheta, lg, ug, options, g = Theta_setup(Prev_Run, num_sec, ng)

        thetaopt, fopt, info = minimize(wing_optimizer, theta0, ng, ltheta, utheta, lg, ug, options)

        println("Optimized twist values: ", thetaopt*180/pi)

        if Include_Airframe == 1
            
            span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, dh, dv = AirFrameVariables(scale_factor, num_sec, Chord_Length)
            D, weight, rho, Vinf, Sref, CL, chords, fs = Airframe_Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

        else

            span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Variables(scale_factor, num_sec, Chord_Length)
            D, weight, rho, Vinf, Sref, CL, chords, fs = Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

        end

        xle_opt = xle
        theta = thetaopt

        #Reference Parameters
        Sref, cref, rref = GetRef(c, chords, theta, span, num_sec, yle)

        ref = Reference(Sref, cref, span, rref, Vinf)

        if Include_Airframe == 1

            system_opt, grids, surfaces_opt, properties_opt, symmetric = AircraftPaneling(xle_opt, yle, zle, chords, theta, phi, ns, nc,
            spacing_s, spacing_c, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h, mirror_h, fc_h, spacing_s_h, spacing_c_h, xle_v, 
            yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v, mirror_v, fc_v, spacing_s_v, spacing_c_v, dh, dv)
            
            properties_opt = get_surface_properties(system_opt)
            
            r, c = lifting_line_geometry(grids, 0.25)
            
            cf, cm = lifting_line_coefficients(system_opt, r, c; frame=Wind())
            
            write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

            yle_2, y2, Lift_prime, elliptical_distribution = CreateLiftCoefficients(cf, span, rho, Vinf, chords, num_sec)

        else
            # reconstruct surface
            grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chords, thetaopt, phi, ns, nc;
            spacing_s=spacing_s, spacing_c=spacing_c)

            surfaces_opt, system_opt, properties_opt = VTK_setup(surface_opt, Vinf, alpha_angle, beta, Omega, ref, symmetric)

            write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

            frame=Wind()

            cf = LocalizedCoefficients(system_opt, grid_opt, frame)
            
            yle_2, y2, Lift_prime, elliptical_distribution = GetLiftCurves(cf, chords, rho, Vinf, num_sec, span)

        end

        PlotLiftDistr(yle_2, y2, Lift_prime, elliptical_distribution)   

return thetaopt
end

function optimizechord(;chordopt)

        num_sec = 12
        scale_factor = 1

        if chordopt != 0
        Prev_Run = 1
        else
        Prev_Run = 0
        end
        
        Include_Airframe = 0
        ng = 2*num_sec+1

        #Creating the optimization problem
        function wing_optimizer(g, chords)

            chordopt=chords

            if Include_Airframe == 1
            
                D, weight, rho, Vinf, Sref, CL, theta, fs = Airframe_Drag_Calculation_Chord(chordopt, num_sec, scale_factor)

            else

                D, weight, rho, Vinf, Sref, CL, theta, fs = Drag_Calculation_Chord(chordopt, num_sec, scale_factor)

            end

            g[1]=weight-.5*rho*Vinf^2*Sref*CL

            # Calculate chord differences
            for i in 1:num_sec
                g[i+1] = chordopt[i+1] - chordopt[i]#+.25/num_sec
            end

        return D
        end

        chord0, lchord, uchord, lg, ug, options, g = chord_setup(Prev_Run, num_sec, ng)

        chordopt, fopt, info = minimize(wing_optimizer, chord0, ng, lchord, uchord, lg, ug, options)

        println("Optimized chord values: ", chordopt)

        if Include_Airframe == 1
            
            span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, dh, dv = AirFrameChordVariables(scale_factor, num_sec, chordopt)
            D, weight, rho, Vinf, Sref, CL, chords, fs = Airframe_Drag_Calculation_Chord(chordopt, num_sec, scale_factor)

        else

            span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Chord_Variables(scale_factor, num_sec, chordopt)
            D, weight, rho, Vinf, Sref, CL, chords, fs = Drag_Calculation_Chord(chordopt, num_sec, scale_factor)

        end

        #Reference Parameters
        Sref, cref, rref = GetRef(c, chords, theta, span, num_sec, yle)

        ref = Reference(Sref, cref, span, rref, Vinf)

        xle_opt = xle

        # Plot the chords
        plot_chords(xle_opt, yle, chordopt)

        if Include_Airframe == 1

            system_opt, grids, surfaces_opt, properties_opt, symmetric = AircraftPaneling(xle_opt, yle, zle, chords, theta, phi, ns, nc,
            spacing_s, spacing_c, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h, mirror_h, fc_h, spacing_s_h, spacing_c_h, xle_v, 
            yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v, mirror_v, fc_v, spacing_s_v, spacing_c_v, dh, dv)
            
            properties_opt = get_surface_properties(system_opt)
            
            r, c = lifting_line_geometry(grids, 0.25)
            
            cf, cm = lifting_line_coefficients(system_opt, r, c; frame=Wind())
            
            write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

            yle_2, y2, Lift_prime, elliptical_distribution = CreateLiftCoefficients(cf, span, rho, Vinf, chords, num_sec)

        else
            # reconstruct surface
            grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chords, thetaopt, phi, ns, nc;
            spacing_s=spacing_s, spacing_c=spacing_c)

            surfaces_opt, system_opt, properties_opt = VTK_setup(surface_opt, Vinf, alpha_angle, beta, Omega, ref, symmetric)

            write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=symmetric)

            frame=Wind()

            cf = LocalizedCoefficients(system_opt, grid_opt, frame)
            
            yle_2, y2, Lift_prime, elliptical_distribution = GetLiftCurves(cf, chords, rho, Vinf, num_sec, span)

        end

        PlotLiftDistr(yle_2, y2, Lift_prime, elliptical_distribution)   

return chordopt
end