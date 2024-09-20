using VortexLattice
using Plots

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false

# horizontal stabilizer
xle_h = [0.0, 0.14]
yle_h = [0.0, 1.25]
zle_h = [0.0, 0.0]
chord_h = [0.7, 0.42]
theta_h = [0.0, 0.0]
phi_h = [0.0, 0.0]
fc_h = fill((xc) -> 0, 2) # camberline function for each section
ns_h = 6
nc_h = 3
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
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false

Sref = 9.0
cref = 0.9
bref = 10.0
rref = [0.5, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

my_alpha = 5.0*pi/180  # Use a different variable name
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, my_alpha, beta, Omega)

symmetric = [true, true, false]

# generate surface panels for wing
wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

# generate surface panels for horizontal tail
hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
VortexLattice.translate!(hgrid, [4.0, 0.0, 0.0])
VortexLattice.translate!(htail, [4.0, 0.0, 0.0])

# generate surface panels for vertical tail
vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
VortexLattice.translate!(vgrid, [4.0, 0.0, 0.0])
VortexLattice.translate!(vtail, [4.0, 0.0, 0.0])

# now set normal vectors manually
ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

# overwrite normal vector for each wing panel
for i = 1:length(wing)
    wing[i] = set_normal(wing[i], ncp)
end

grids = [wgrid, hgrid, vgrid]
surfaces = [wing, htail, vtail]
surface_id = [1, 2, 3]

system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

CF, CM = body_forces(system; frame=Stability())

CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

properties = get_surface_properties(system)

write_vtk("wing-tail", surfaces, properties; symmetric)