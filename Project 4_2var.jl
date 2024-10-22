using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots

global num_sec = 12
global sec_2 = Int(.5*num_sec)
global scale_factor = 1
global alpha_max = 25

#Creating the optimization problem
function wing_optimizer(g, x)

    # Define inputs of function
    span = 8.0 #one wing or the whole span       
    rho = 1.225
    weight = 1.7*scale_factor
    Vinf = 1.0

    # geometry (right half of the wing)
    yle = [i * (span / (num_sec)) for i in 0:(num_sec)]
    zle = zeros(num_sec+1)
    theta = zeros(num_sec+1)
    phi = zeros(num_sec+1)
    chords = zeros(num_sec+1)

    for i in 1:num_sec+1
        chords[i] = x[i]
    end

    for i in num_sec+2:2*num_sec+2
        theta[i-(num_sec+1)] = x[i]
    end

    fc = zeros(num_sec+1)
    xle = zeros(num_sec+1)

    for i in 1:num_sec
        xle[i+1] = (chords[1]/4 - chords[i+1]/4)
    end


    # discretization parameters
    ns = num_sec
    nc = num_sec
    spacing_s = Uniform()
    spacing_c = Uniform()

    Sref= 0.0

    # Reference Area Calculation
    for i in 1:num_sec
        S = ((chords[i] + chords[i+1]) / 2) * (yle[i+1] - yle[i])
        Sref=S+Sref
    end

    cref= Sref/span

    # reference parameters
    rref = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, span, rref, Vinf)


    # freestream parameters
    alpha_angle = 5*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha_angle, beta, Omega)

    # construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces = [surface]

    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = true

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

    g[1]=weight-.5*rho*Vinf^2*Sref*CL
    g[2]=x[1]-x[2]-.02
+i]-x[num_sec+1+i]
        println("Negative")
        else
        g[i+2*num_sec+2]=x[num_sec+1+i]-x[num_sec+1+i+1]
        end
    end

    return D
end

# Initialize vectors based on num_sec
x0 = zeros(2*num_sec+2)  # starting point

# Initialize vectors based on num_sec
for i in 1:num_sec+1
x0[i] = 1  # starting point
end

if num_sec >= length(chord_opt)
for i in 1:length(chord_opt)
x0[i]=xopt[i]
end
else
    for i in 1:num_sec+1
    x0[i]=xopt[i]
    end
end

if num_sec >= length(chord_opt)
    for i in length(chord_opt)+1:2*length(chord_opt)+1
    x0[i]=xopt[i]
    end
    else
        for i in num_sec+2:2*num_sec+2
        x0[i]=xopt[i]
    end
end

lx = fill(-45.0*pi/180, 2*num_sec+2)  # lower bounds on x
ux = fill(45.0*pi/180, 2*num_sec+2)  # upper bounds on x

# Initialize vectors based on num_sec
for i in 1:num_sec+1
    lx[i] = 0.01  # starting point
end

# Initialize vectors based on num_sec
for i in 1:num_sec+1
    ux[i] = 5.0  # starting point
end

ng = 2 + 3*num_sec  # number of constraints
lg = -Inf*one(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g
g = zeros(ng)

# ----- set some options ------
ip_options = Dict(
    "max_iter" => 150,
    "tol" => .9
)
solver = IPOPT(ip_options)
options = Options(;solver, derivatives=ForwardFD())

xopt, fopt, info = minimize(wing_optimizer, x0, ng, lx, ux, lg, ug, options)

span = 8.0 #one wing or the whole span       
rho = 1.225
Vinf = 1.0

chord_opt=zeros(num_sec+1)

for i in  1:num_sec+1
chord_opt[i]=xopt[i]
end

xle_opt = zeros(num_sec+1)

for i in 1:num_sec
    xle_opt[i+1] = (chord_opt[1]/4 - chord_opt[i+1]/4)
end

for i in num_sec+2:2*num_sec+2
    theta[i-(num_sec+1)] = xopt[i]
end

# discretization parameters
ns = num_sec
nc = num_sec

symmetric = true

# geometry (right half of the wing)
yle = [i * (span / (num_sec)) for i in 0:(num_sec)]
zle = zeros(num_sec+1)
theta = zeros(num_sec+1)
phi = zeros(num_sec+1)

spacing_s = Uniform()
spacing_c = Uniform()

Sref= 0.0

# Reference Area Calculation
for i in 1:num_sec
    global Sref
    S = ((chord_opt[i] + chord_opt[i+1]) / 2) * (yle[i+1] - yle[i])
    Sref=S+Sref
end

cref= Sref/span

# reference parameters
rref = [0.50, 0.0, 0.0]
ref = Reference(Sref, cref, span, rref, Vinf)

alpha_opt=5

# freestream parameters
alpha_angle = alpha_opt*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha_angle, beta, Omega)

# reconstruct surface
grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chord_opt, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces_opt = [surface_opt]

system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

properties_opt = get_surface_properties(system_opt)

grids=[grid_opt]

r, c = lifting_line_geometry(grids, 0.25)

cf, cm = lifting_line_coefficients(system_opt, r, c; frame=Wind())

write_vtk("optimized-symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=true)

println("Optimized leading edge values: ", xle_opt)
println("Optimized chord values: ", chord_opt)
# println("Optimized alpha value: ", alpha_opt)
println("Optimized twist values: ", thetaopt*180/pi)

# Plotting function
function plot_chords(xle_opt, yle, chords)
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

# Plot the chords
plot_chords(xle_opt, yle, chord_opt)

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
savefig("Lift_Distribution_along_the_Span_Twist_Optimization.pdf")