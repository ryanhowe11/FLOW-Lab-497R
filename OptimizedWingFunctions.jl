    #Setting up all the problem variables
    function ProblemSetup(num_sec, scale_factor, chord_opt)
        # Define inputs of function
     span = 8.0 #one wing or the whole span       
     rho = 1.225
     Vinf = 1.0
     weight = 1.7*scale_factor

     T = eltype(chord_opt)

	chords = zeros(T, num_sec+1)

   	 for i in 1:num_sec+1
     chords[i] = chord_opt[i]
     end

     Vinf=5*chord_opt[num_sec+2]

     # geometry (right half of the wing)
     yle = zeros(T, num_sec+1)
     yle = [i * (span / (num_sec)) for i in 0:(num_sec)]
     zle = zeros(T, num_sec+1)
     theta = zeros(T, num_sec+1)
     phi = zeros(T, num_sec+1)
     fc = zeros(T, num_sec+1)
     xle = zeros(T, num_sec+1)

         # discretization parameters
     ns = num_sec
     nc = num_sec
     spacing_s = Uniform()
     spacing_c = Uniform()
    
 return chords, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c
 end

 #Find the leading edge values assuming zero sweep
 function XLEcalc(xle, chords, num_sec)

     for i in 1:num_sec
         xle[i+1] = (chords[1]/4 - chords[i+1]/4)
     end

     return xle
 end

 #Define your references such as area and position
 function ReferenceCalculation(num_sec, yle, span, Vinf, c)

     Sref= 0.0

     # Reference Area Calculation
     for i in 1:num_sec
         S = ((c[i] + c[i+1]) / 2) * (yle[i+1] - yle[i])
         Sref=S+Sref
     end

     cref= Sref/span

     # reference parameters
     rref = [0.50, 0.0, 0.0]
     ref = Reference(Sref, cref, span, rref, Vinf)

 return Sref, ref
 end

 #gives you your fs variable to use in Vortex Lattice
 function FreestreamParams(Vinf)

     # freestream parameters
     alpha_angle = 5.0*pi/180
     beta = 0.0
     Omega = [0.0; 0.0; 0.0]
     fs = Freestream(Vinf, alpha_angle, beta, Omega)

 return fs
 end

 #Finds the coefficients of lift and drag for whole system
 function CalculateSurfaceForcesCoeff(surface, ref, fs, symmetric)

     # create vector containing all surfaces
     surfaces = [surface]

     # perform steady state analysis
     system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

     # retrieve near-field forces
     CF, CM = body_forces(system; frame=Wind())

     CD, CY, CL = CF
     Cl, Cm, Cn = CM

 return CL, CD
 end

 #Setting Up Optimization Problem Variables
 function OptimizationSetup(num_sec)
     # Initialize vectors based on num_sec
     c0 = ones(num_sec+2)  # starting point
 
     # if num_sec >= length(chord_opt)
     # for i in 1:length(chord_opt)
     # c0[i]=chord_opt[i]
     # end
     # else
     #     for i in 1:num_sec
     #     c0[i]=chord_opt[i]
     #     end
     # end
 
     lc = fill(0.01, num_sec+2)  # lower bounds on x
     uc = fill(5.0, num_sec+2)  # upper bounds on x
     ng = 2 + num_sec  # number of constraints
     lg = -Inf*one(ng)  # lower bounds on g
     ug = zeros(ng)  # upper bounds on g
     g = zeros(ng)
 
     # ----- set some options ------
     ip_options = Dict(
         "max_iter" => 100,
         "tol" => 1e-3
     )
     solver = IPOPT(ip_options)
     options = Options(;solver, derivatives=ForwardAD())
     return c0, ng, lc, uc, lg, ug, options, g
 end
 
 # Plotting function
 function plot_chords(xle_opt, yle, chords, num_sec)
     plot()
     for i in 1:num_sec
         x_start = xle_opt[i]
         y_start = yle[i]
         x_end = x_start + chords[i]
         y_end = y_start
         plot!([x_start, x_end], [y_start, y_end], label="Chord $i", legend=false)
     end
     xlabel!("yle")
     ylabel!("xle_opt / Chord length")
     title!("Chords at Each Section")
 end

 function RunOptimizer(num_sec, scale_factor)

    c0, ng, lc, uc, lg, ug, options, g = OptimizationSetup(num_sec)

    #Creating the optimization problem
    function wing_optimizer(g, c)

        chords, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, scale_factor, c)

        xle = XLEcalc(xle, chords, num_sec)

        Sref, ref = ReferenceCalculation(num_sec, yle, span, Vinf, c)

        fs = FreestreamParams(Vinf)

        # construct surface
        grid, surface = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

            # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
        symmetric = true

        CL, CD = CalculateSurfaceForcesCoeff(surface, ref, fs, symmetric)

        D=.5*rho*Vinf^2*Sref*CD
        L = .5*rho*Vinf^2*Sref*CL

        g[1]=weight-L
        g[2]=c[1]-c[2]-.02

        # Calculate chord differences
        for i in 1:num_sec
            g[i+2] = c[i+1] - c[i]
        end

        O=D/(L*Vinf)

    return O
    end

    xopt, fopt, info = minimize(wing_optimizer, c0, ng, lc, uc, lg, ug, options)

return xopt, fopt, info
end