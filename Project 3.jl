using Plots
using VortexLattice
using QuadGK
using Interpolations
using Trapz
using SNOW

# Define inputs of function
span = 4 #one wing or the whole span       
rho = 1.225
weight = 1.7
Vinf = 1.0

#Creating the optimization problem
function optimize_wing()
    model = Model(Ipopt.Optimizer)

    #Define variables
    @variable(model, x[1:7], start = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    @variable(model, c[1:7], start = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4])

    #Define constraints for the Optimizer
    for i in 1:6
        @constraint(model, x[i+1] > x[i])
        @constraint(model, c[i] > c[i+1])
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

# perform far-field analysis
CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

dCF, dCM = stability_derivatives(system)

CDa, CYa, CLa = dCF.alpha
Cla, Cma, Cna = dCM.alpha
CDb, CYb, CLb = dCF.beta
Clb, Cmb, Cnb = dCM.beta
CDp, CYp, CLp = dCF.p
Clp, Cmp, Cnp = dCM.p
CDq, CYq, CLq = dCF.q
Clq, Cmq, Cnq = dCM.q
CDr, CYr, CLr = dCF.r
Clr, Cmr, Cnr = dCM.r

properties = get_surface_properties(system)

write_vtk("symmetric-planar-wing", surfaces, properties; symmetric)

grids=[grid]

r, c = lifting_line_geometry(grids, 0.25)

cf, cm = lifting_line_coefficients(system, r, c; frame=Wind())

A = 0.0  # Initialize the area

# Ensure A is treated as a local variable within the loop
for i = 1:6
    global A
    A += ((chord[i] + chord[i+1]) / 2) * (yle[i+1] - yle[i])
end

L=.5*rho*Vinf^2*A*CL
D=.5*rho*Vinf^2*A*CD