function Variables(scale_factor)

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

function ReferenceAreaTheta(c, chords, theta)

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

function GetRef(c, chords, theta, span)

    Sref = ReferenceAreaTheta(c, chords, theta)

    cref= Sref/span

    # reference parameters
    rref = [0.50, 0.0, 0.0]

    return Sref, cref, rref
end

function Drag_Calculation_Theta(thetaopt)

    span, rho, weight, Vinf, yle, zle, chords, theta, phi, fc, xle, c, ns, nc, spacing_s, spacing_c, alpha_angle, beta, Omega, symmetric = Variables()

    theta=thetaopt

    Sref, cref, rref = GetRef(c, chords, theta, span)
    
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

return D, weight, rho, Vinf, Sref, CL, chords
end

function setup(Prev_Run)

    # Initialize vectors based on num_sec
    theta0 = fill(5*pi/180, num_sec+1)  # starting point

    if Prev_Run == 1
        if num_sec >= length(thetaopt)
        for i in 1:length(thetaopt)
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
    ng = 1 #+ num_sec  # number of constraints
    lg = -Inf*ones(ng)  # lower bounds on g
    ug = zeros(ng)  # upper bounds on g
    g = zeros(ng)

    # ----- set some options ------
    ip_options = Dict(
        "max_iter" => 100,
        "tol" => 1e-6
    )
    solver = IPOPT(ip_options)
    options = Options(;solver, derivatives=ForwardFD())

return theta0, ng, ltheta, utheta, lg, ug, options, g
end

function LocalizedCoefficients(system_opt, grid_opt, frame)

    grids=[grid_opt]
    
    r, c = lifting_line_geometry(grids, 0.25)
    
    cf, cm = lifting_line_coefficients(system_opt, r, c; frame=frame)

return cf
end

function XLE_calc(chords, num_sec)

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