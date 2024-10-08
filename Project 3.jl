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
yle = (0, 2/3, 4/3, 2, 8/3, 10/3, 4)
zle = (0, 0, 0, 0, 0, 0, 0)
theta = fill(0, 7)
phi = zeros(7)
chords=c
fc = fill((xc) -> 0, 7)                     # camberline function for each section

xle=zeros(7)

for i in 1:6
xle[i+1]=c[i]/4-c[i+1]/4
end


# discretization parameters
ns = 7
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

for i in 1:6
g[1+i]=x[i]-x[i+1]
g[8+i]=c[i+1]-c[i]
end

return D
end

c0 = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]  # starting point
lc = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]  # lower bounds on x
uc = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]  # upper bounds on x
ng = 13  # number of constraints
lg = -Inf*one(ng)  # lower bounds on g
ug = zeros(ng)  # upper bounds on g

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