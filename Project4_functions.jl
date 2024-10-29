function Variables(scale_factor, num_sec, Chord_Length)

    # Define inputs of function
    span = 8.0 #one wing or the whole span       
    rho = 1.225
    weight = 1.7*scale_factor
    Vinf = 1.0

    # geometry (right half of the wing)
    yle = [i * (span / (num_sec)) for i in 0:(num_sec)]
    zle = zeros(num_sec+1)
    chords = ones(num_sec+1)*Chord_Length
    theta = zeros(num_sec+1)
    phi = zeros(num_sec+1)
    fc = zeros(num_sec+1)
    xle = zeros(num_sec+1)
    c=zeros(num_sec+1)

    # discretization parameters
    ns = num_sec+1
    nc = num_sec
    spacing_s = Uniform()
    spacing_c = Uniform()

    # freestream parameters
    alpha_angle = 5*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

return span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric
end

function AirFrameVariables(scale_factor, num_sec, Chord_Length)

    # Define inputs of function
    span = 8.0 #one wing or the whole span       
    rho = 1.225
    weight = 1.7*scale_factor
    Vinf = 1.0

    dh = 4
    dv = 4
    
    # horizontal stabilizer
    xle_h = [0.0, 0.14]
    yle_h = [0.0, 1.25]
    zle_h = [0.0, 0.0]
    chord_h = [0.7, 0.42]
    theta_h = [0.0, 0.0]
    phi_h = [0.0, 0.0]
    fc_h = fill((xc) -> 0, 2) # camberline function for each section
    ns_h = num_sec
    nc_h = 1
    spacing_s_h = Uniform()
    spacing_c_h = Uniform()
    mirror_h = false
    
    # vertical stabilizer
    xle_v = [0.0, 0.14]
    yle_v = [0.0, 0.0]
    zle_v = [0.0, 1.0]
    chord_v = [0.7, 0.42]
    theta_v = [0.0, 0.0]
    phi_v = [0.0, 0.0]
    fc_v = fill((xc) -> 0, 2) # camberline function for each section
    ns_v = num_sec
    nc_v = 1
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    # geometry (right half of the wing)
    yle = [i * (span / (num_sec)) for i in 0:(num_sec)]
    zle = zeros(num_sec+1)
    chords = ones(num_sec+1)*Chord_Length
    theta = zeros(num_sec+1)
    phi = zeros(num_sec+1)
    fc = zeros(num_sec+1)
    xle = zeros(num_sec+1)
    c=zeros(num_sec+1)

    # discretization parameters
    ns = num_sec+1
    nc = num_sec
    spacing_s = Uniform()
    spacing_c = Uniform()

    # freestream parameters
    alpha_angle = 5*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

return span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, dh, dv
end

function ReferenceAreaTheta(c, chords, theta, num_sec, yle)

    Sref= 0.0

    for i in 1:num_sec
        c[i]=chords[i]*cos(theta[i])
    end

    # Reference Area Calculation
    for i in 1:num_sec
        S = ((c[i] + c[i+1]) / 2) * (yle[i+1] - yle[i])
        Sref=S+Sref
    end

return Sref
end

function GetRef(c, chords, theta, span, num_sec, yle)

    Sref = ReferenceAreaTheta(c, chords, theta, num_sec, yle)

    cref= Sref/span

    # reference parameters
    rref = [0.50, 0.0, 0.0]

    return Sref, cref, rref
end

function Airframe_Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

    span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v, dh, dv = AirFrameVariables(scale_factor, num_sec, Chord_Length)

    theta=thetaopt

    Sref, cref, rref = GetRef(c, chords, theta, span, num_sec, yle)
    
    ref = Reference(Sref, cref, span, rref, Vinf)

    fs = Freestream(Vinf, alpha_angle, beta, Omega)

    # construct surface
    wgrid, wing = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
            spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    VortexLattice.translate!(hgrid, [dh, 0.0, 0.0])
    VortexLattice.translate!(htail, [dh, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    VortexLattice.translate!(vgrid, [dv, 0.0, 0.0])
    VortexLattice.translate!(vtail, [dv, 0.0, 0.0])

    # create vector containing all surfaces
    grids = [wgrid, hgrid, vgrid]
    surfaces = [wing, htail, vtail]
    surface_id = [1, 2, 3]

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    println(CF)

    # perform far-field analysis
    CDiff = far_field_drag(system)

    dCF, dCM = stability_derivatives(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    properties = get_surface_properties(system)

    #write_vtk("intermediate-symmetric-planar-wing", surfaces, properties; symmetric)

    D=.5*rho*Vinf^2*Sref*CD

return D, weight, rho, Vinf, Sref, CL, chords, fs
end

function Drag_Calculation_Theta(thetaopt, num_sec, scale_factor, Chord_Length)

    span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Variables(scale_factor, num_sec, Chord_Length)

    theta=thetaopt

    Sref, cref, rref = GetRef(c, chords, theta, span, num_sec, yle)
    
    ref = Reference(Sref, cref, span, rref, Vinf)

    fs = Freestream(Vinf, alpha_angle, beta, Omega)

    # construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces = [surface]

    # perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    println(CF)

    # perform far-field analysis
    CDiff = far_field_drag(system)

    dCF, dCM = stability_derivatives(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    properties = get_surface_properties(system)

    #write_vtk("intermediate-symmetric-planar-wing", surfaces, properties; symmetric)

    D=.5*rho*Vinf^2*Sref*CD

return D, weight, rho, Vinf, Sref, CL, chords, fs
end

function setup(Prev_Run, num_sec, ng)

    # Initialize vectors based on num_sec
    theta0 = fill(5*pi/180, num_sec+1)  # starting point

    if Prev_Run == 1
        if num_sec >= length(thetaopt)
        for i = 1:length(thetaopt)
        theta0[i]=thetaopt[i]
        end
        else
            for i in 1:num_sec
            theta0[i]=thetaopt[i]
            end
        end
    end

    ltheta = fill(-45.0*pi/180, num_sec+1)  # lower bounds on x
    utheta = fill(45.0*pi/180, num_sec+1)  # upper bounds on x
    lg = -Inf*ones(ng)  # lower bounds on g
    ug = zeros(ng)  # upper bounds on g
    g = zeros(ng)

    # ----- set some options ------
    ip_options = Dict(
        "max_iter" => 100,
        "tol" => 1e-3
    )
    solver = IPOPT(ip_options)
    options = Options(;solver, derivatives=ForwardFD())

return theta0, ltheta, utheta, lg, ug, options, g
end

function LocalizedCoefficients(system_opt, grid_opt, frame)

    grids=[grid_opt]
    
    r, c = lifting_line_geometry(grids, 0.25)
    
    cf, cm = lifting_line_coefficients(system_opt, r, c; frame=frame)

return cf
end

function XLE_calc(chords, num_sec)

    xle_opt = zeros(num_sec+1)

    for i in 1:num_sec
        xle_opt[i] = ((1-chords[i])/2)
    end

return xle_opt
end

function GetLiftCurves(cf, chords, rho, Vinf, num_sec, span)

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

    Lift_prime = .5*rho*Vinf^2 .* chords .* lifting_coefficients[:, 1]

    yle_2 = [i * (span / (num_sec)) for i in 0:(num_sec)]

    println(length(yle_2), length(Lift_prime))

    area_prime = trapz(yle_2, Lift_prime)

    println(area_prime)

    bprime = span
    aprime = 4*area_prime/(bprime*pi)

    # Calculate the ideal elliptic lift distribution
    θ = range(0, π/2, length=100)
    x = bprime * cos.(θ)
    y = aprime * sin.(θ)
    Cl_max = maximum(Lift_prime)
    elliptical_distribution = Cl_max * sqrt.(1 .- (y2 ./ span).^2)

return yle_2, y2, Lift_prime, elliptical_distribution
end

function PlotLiftDistr(yle_2, y2, Lift_prime, elliptical_distribution)
    # Plot the lift distribution
    plot(yle_2, Lift_prime, label="Optimized Lift Distribution", xlabel="Spanwise Location (y)", ylabel="Lift Coefficient (Cl)")
    plot!(y2, elliptical_distribution, label="Elliptical Lift Distribution", linestyle=:dash)

    # Save the lift distribution plot as a PDF
    savefig("Lift_Distribution_along_the_Span_Twist_Optimization.pdf")
end

function VTK_setup(surface_opt, Vinf, alpha_angle, beta, Omega, ref, symmetric)

    # create vector containing all surfaces
    surfaces_opt = [surface_opt]

    fs = Freestream(Vinf, alpha_angle, beta, Omega)

    system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

    properties_opt = get_surface_properties(system_opt)

return surfaces_opt, system_opt, properties_opt
end

# Define the avl_normal_vector function
function avl_normal_vector(vector, angle)
    # Example calculation
    return normalize(vector) * cos(angle)
end

function AircraftPaneling(xle_opt, yle, zle, chords, theta, phi, ns, nc,
    spacing_s, spacing_c, xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h, mirror_h, fc_h, spacing_s_h, spacing_c_h, xle_v, 
    yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v, mirror_v, fc_v, spacing_s_v, spacing_c_v, dh, dv)

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

        return system_opt, grids, surfaces_opt, properties_opt, symmetric
end

function CreateLiftCoefficients(cf, span, rho, Vinf, chords, num_sec)

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

    return yle_2, y2, Lift_prime, elliptical_distribution
end