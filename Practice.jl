
using Plots
using VortexLattice
using QuadGK
using Interpolations
using Trapz
using SNOW
using Ipopt


# Define inputs of function
span = 4 #one wing or the whole span       
rho = 1.225
weight = 1.7
Vinf = 1.0

#Creating the optimization problem
function wing_optimizer()
    model = Model(Ipopt.Optimizer)

    #Define variables
    @variable(model, x[1:7])
    @variable(model, c[1:7])

    for i in 1:6
        set_start_value(x[i], 0)
        set_start_value(c[i], 1)
    end

    #=
    @variable(model, x[1] == 0.0)
    @variable(model, x[2] == 0.0)
    @variable(model, x[3] == 0.0)
    @variable(model, x[4] == 0.0)
    @variable(model, x[5] == 0.0)
    @variable(model, x[6] == 0.0)
    @variable(model, x[7] == 0.0)

    @variable(model, c[1] == 4.0)
    @variable(model, c[2] == 4.0)
    @variable(model, c[3] == 4.0)
    @variable(model, c[4] == 4.0)
    @variable(model, c[5] == 4.0)
    @variable(model, c[6] == 4.0)
    @variable(model, c[7] == 4.0)
=#

    #Define constraints for the Optimizer
    for i in 1:6
        @constraint(model, x[i+1] >= x[i])
        @constraint(model, c[i] >= c[i+1])
    end

# geometry (right half of the wing)
xle = x
yle = (0, 2/3, 4/3, 2, 8/3, 10/3, 4)
zle = (0, 0, 0, 0, 0, 0 ,0)
chord = c
theta = fill(0, 7)
phi = zeros(7)
fc = fill((xc) -> 0, 7)                     # camberline function for each section

# discretization parameters
ns = 7
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

Sref= 0.0

# Reference Area Calculation
for i in 1:6
    global Sref
    Sref += ((chord[i] + chord[i+1]) / 2) * (yle[i+1] - yle[i])
end

# reference parameters
rref = [0.50, 0.0, 0.0]
ref = Reference(Sref, x[1], span, rref, Vinf)


# freestream parameters
alpha_angle = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha_angle, beta, Omega)

# construct surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces = [surface]

# we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
symmetric = true

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

CD, CY, CL = CF
Cl, Cm, Cn = CM

properties = get_surface_properties(system)

write_vtk("optimized-symmetric-planar-wing", surfaces, properties; symmetric)

L=.5*rho*Vinf^2*Sref*CL
D=.5*rho*Vinf^2*Sref*CD

#Set forth objective to minimize drag
@objective(model, Min, D)

#Lift constraint
@constraint(model, L >= weight)

# Solve the optimization problem
optimize!(model)

# Extract optimized values
x_opt = value.(x)
c_opt = value.(c)

#update the geometry with optimized values
xle_opt = x_opt
chord_opt = c_opt

# reconstruct surface
grid_opt, surface_opt = wing_to_surface_panels(xle_opt, yle, zle, chord_opt, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces_opt = [surface_opt]

system = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

properties_opt = get_surface_properties(system)

write_vtk("symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=true)

return x_opt, c_opt
end

x_opt, c_opt = wing_optimizer()
println("Optimized leading edge values: ", x_opt)
println("Optimized chord values: ", c_opt)