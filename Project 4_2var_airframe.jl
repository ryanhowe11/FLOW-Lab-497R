using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff
using Plots
using LinearAlgebra

global num_sec = 12
global sec_2 = Int(.5*num_sec)
global scale_factor = 1
global alpha_max = 25
global Chord_Length = 1
#Creating the optimization problem
function wing_optimizer(g, x)

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
    theta = zeros(num_sec+1)
    phi = zeros(num_sec+1)
    chords = zeros(num_sec+1)

    for i in 1:num_sec+1
    chords[i] = x[i]
    end

    fc = zeros(num_sec+1)
    xle = zeros(num_sec+1)

    for i in 1:num_sec
        xle[i+1] = (chords[1]/4 - chords[i+1]/4)
    end


    # discretization parameters
    ns = num_sec
    nc = 1
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
    alpha_angle = x[num_sec+2]*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
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

# now set normal vectors manually
ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

# overwrite normal vector for each wing panel
for i = 1:length(wing)
    wing[i] = set_normal(wing[i], ncp)
end

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
    L=.5*rho*Vinf^2*Sref*CL

    Einv=(D/L)

    g[1]=weight-L
    g[2]=-1

    # Calculate xle differences
    for i in 1:num_sec
        g[i+2] = xle[i] - xle[i+1]
    end

    # Calculate chord differences
    for i in 1:num_sec
        g[i+num_sec+2] = x[i+1] - x[i]
    end


    for i in 1:sec_2
    g[i+2+2*num_sec]=x[i]-x[i+1]-.5
    end


    return Einv
end

# function to construct a normal vector the way AVL does
#  - `ds` is a line representing the leading edge
#  - `theta` is the incidence angle, taken as a rotation (+ by RH rule) about
#        the surface's spanwise axis projected onto the Y-Z plane.
function avl_normal_vector(ds, theta)

    st, ct = sincos(theta)

    # bound vortex vector
    bhat = ds/norm(ds)

    # chordwise strip normal vector
    shat = [0, -ds[3], ds[2]]/sqrt(ds[2]^2+ds[3]^2)

    # camberline vector
    chat = [ct, -st*shat[2], -st*shat[3]]

    # normal vector perpindicular to camberline and bound vortex for entire chordwise strip
    ncp = cross(chat, ds)
    return ncp / norm(ncp) # normal vector used by AVL
end

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

# Initialize vectors based on num_sec
x0 = ones(num_sec+2)  # starting point

# if num_sec >= length(chord_opt)
# for i in 1:length(chord_opt)
# c0[i]=chord_opt[i]
# end
# else
#     for i in 1:num_sec
#     c0[i]=chord_opt[i]
#     end
# end


lx = fill(0.01, num_sec+2)  # lower bounds on x
ux = fill(5.0, num_sec+2)  # upper bounds on x
ux[num_sec+2]= alpha_max
ng = 2 + sec_2 + 2*num_sec  # number of constraints
lg = -Inf*one(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g
g = zeros(ng)

# ----- set some options ------
ip_options = Dict(
    "max_iter" => 100,
    "tol" => 1e-3
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

# discretization parameters
ns = num_sec
nc = 1

symmetric = [true, true, false]

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

alpha_opt=xopt[num_sec+2]

# freestream parameters
alpha_angle = alpha_opt*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha_angle, beta, Omega)

# reconstruct surface
wgrid_opt, wing_opt = wing_to_surface_panels(xle_opt, yle, zle, chord_opt, theta, phi, ns, nc;
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
ncp = avl_normal_vector([xle_opt[2]-xle_opt[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

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

println("Optimized leading edge values: ", xle_opt)
println("Optimized chord values: ", chord_opt)
println("Optimized alpha value: ", alpha_opt)

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
savefig("Chords_plot_airframe.pdf")

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