using VortexLattice

# geometry (right half of the wing)
xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 0.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section

# discretization parameters
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
Sref = 30.0
cref = 2.0
bref = 15.0
rref = [0.50, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

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

properties = get_surface_properties(system)

write_vtk("symmetric-planar-wing", surfaces, properties; symmetric)