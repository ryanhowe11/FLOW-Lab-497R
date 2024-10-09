using VortexLattice
using SNOW
using Ipopt
using FiniteDiff
using ForwardDiff

#Creating the optimization problem
function wing_optimizer(g, c)

# Define inputs of function
span = 4.0 #one wing or the whole span       
rho = 1.225
weight = 1.7
Vinf = 1.0

# geometry (right half of the wing)
yle = [0.0, 0.667, 1.333, 2.0, 2.667, 3.333, 4.0]
zle = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
theta = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
phi = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
chords = c
fc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # camberline function for each section
xle = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

for i in 1:6
    xle[i+1] = chords[i]/4 - chords[i+1]/4
end


# discretization parameters
ns = 6
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

Sref= 0.0

# Reference Area Calculation
for i in 1:6
    global Sref
    Sref += ((c[i] + c[i+1]) / 2) * (yle[i+1] - yle[i])
end

x=xle[1]

# reference parameters
rref = [0.50, 0.0, 0.0]
ref = Reference(Sref, x, span, rref, Vinf)


# freestream parameters
alpha_angle = 5.0*pi/180
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

CD, CY, CL = CF
Cl, Cm, Cn = CM

properties = get_surface_properties(system)

write_vtk("optimized-symmetric-planar-wing", surfaces, properties; symmetric)


D=.5*rho*Vinf^2*Sref*CD

g[1]=.5*rho*Vinf^2*Sref*CL-weight
g[2]=xle[1]-xle[2]
g[3]=xle[2]-xle[3]
g[4]=xle[3]-xle[4]
g[5]=xle[4]-xle[5]
g[6]=xle[5]-xle[6]
g[7]=xle[6]-xle[7]
g[8]=c[2]-c[1]
g[9]=c[3]-c[2]
g[10]=c[4]-c[3]
g[11]=c[5]-c[4]
g[12]=c[6]-c[5]
g[12]=c[7]-c[6]

return D
end

c0 = [1.0; .9; .8; .7; .6; .5; .4]  # starting point
lc = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]  # lower bounds on x
uc = [5.0, 5, 5, 5, 5, 5, 5]  # upper bounds on x
ng = 13  # number of constraints
lg = -Inf*one(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g
g = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

# ----- set some options ------
ip_options = Dict(
    "max_iter" => 3,
    "tol" => 1e-6
)
solver = IPOPT(ip_options)
options = Options(;solver)

xopt, fopt, info = minimize(wing_optimizer, c0, ng, lc, uc, lg, ug, options)

#=
# reconstruct surface
grid_opt, surface_opt = wing_to_surface_panels(x_opt, yle, zle, c_opt, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces_opt = [surface_opt]

system_opt = steady_analysis(surfaces_opt, ref, fs; symmetric=symmetric)

properties_opt = get_surface_properties(system_opt)

write_vtk("symmetric-planar-wing", surfaces_opt, properties_opt; symmetric=true)

x_opt, c_opt = wing_optimizer()
println("Optimized leading edge values: ", x_opt)
println("Optimized chord values: ", c_opt)
=#