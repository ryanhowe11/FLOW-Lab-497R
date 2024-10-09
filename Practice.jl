using Plots
using VortexLattice
using QuadGK
using Interpolations
using Trapz
using SNOW
using Ipopt

# Define inputs of function
span = 4 # one wing or the whole span       
rho = 1.225
weight = 1.7
Vinf = 1.0

# Objective function to minimize drag
function objective(xc)
    x, c = xc[1:7], xc[8:14]
    
    # Geometry (right half of the wing)
    xle = x
    yle = [0, 2/3, 4/3, 2, 8/3, 10/3, 4]
    zle = [0, 0, 0, 0, 0, 0, 0]
    chord = c
    theta = fill(0, 7)
    phi = zeros(7)
    fc = fill((xc) -> 0, 7)  # camberline function for each section

    # Discretization parameters
    ns = 7
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()

    # Reference Area Calculation
    Sref = sum(((chord[i] + chord[i+1]) / 2) * (yle[i+1] - yle[i]) for i in 1:6)

    # Reference parameters
    rref = [0.50, 0.0, 0.0]
    ref = Reference(Sref, x, span, rref, Vinf)

    # Freestream parameters
    alpha_angle = 5.0 * pi / 180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha_angle, beta, Omega)

    # Construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
                                           spacing_s=spacing_s, spacing_c=spacing_c)

    # Create vector containing all surfaces
    surfaces = [surface]

    # Perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=true)

    # Retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF

    # Drag calculation
    D = 0.5 * rho * Vinf^2 * Sref * CD

    return D
end

# Constraint function for lift
function constraints(xc)
    x, c = xc[1:7], xc[8:14]
    
    # Geometry (right half of the wing)
    xle = x
    yle = [0, 2/3, 4/3, 2, 8/3, 10/3, 4]
    zle = [0, 0, 0, 0, 0, 0, 0]
    chord = c
    theta = fill(0, 7)
    phi = zeros(7)
    fc = fill((xc) -> 0, 7)  # camberline function for each section

    # Discretization parameters
    ns = 7
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()

    # Reference Area Calculation
    Sref = sum(((chord[i] + chord[i+1]) / 2) * (yle[i+1] - yle[i]) for i in 1:6)

    # Reference parameters
    rref = [0.50, 0.0, 0.0]
    ref = Reference(Sref, x, span, rref, Vinf)

    # Freestream parameters
    alpha_angle = 5.0 * pi / 180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha_angle, beta, Omega)

    # Construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
                                           spacing_s=spacing_s, spacing_c=spacing_c)

    # Create vector containing all surfaces
    surfaces = [surface]

    # Perform steady state analysis
    system = steady_analysis(surfaces, ref, fs; symmetric=true)

    # Retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF

    # Lift calculation
    L = 0.5 * rho * Vinf^2 * Sref * CL

    return [L - weight]
end

# Creating the optimization problem
function wing_optimizer()
    # Initial guess
    x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    # Lower and upper bounds
    lb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ub = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]

    # Define the optimization problem
    prob = SNOW.OptimizationProblem(objective, constraints, x0, lb, ub)

    # Solve the optimization problem using Ipopt
    result = SNOW.optimize(prob, Ipopt.Optimizer)

    # Extract optimized values
    x_opt = result.minimizer[1:7]
    c_opt = result.minimizer[8:14]

    return x_opt, c_opt
end

x_opt, c_opt = wing_optimizer()
println("Optimized leading edge values: ", x_opt)
println("Optimized chord values: ", c_opt)