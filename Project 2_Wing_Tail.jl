using VortexLattice
using Plots

dh = 4
dv = 4

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

symmetric = [true, true, false]

# generate surface panels for wing
wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

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

grids = [wgrid, hgrid, vgrid]
surfaces = [wing, htail, vtail]
surface_id = [1, 2, 3]

# Function to run the simulation and return coefficients
#function get_coefficients(alpha)
fs = Freestream(Vinf, my_alpha, beta, Omega)
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)
CF, CM = body_forces(system; frame=Stability())
CDiff = far_field_drag(system)
CD, CY, CL = CF
Cl, Cm, Cn = CM
#return CD, CY, CL, Cl, Cm, Cn
#end

#=
# Define a small variation in angle of attack (1 degree in radians)
delta_alpha = 1.0 * pi / 180

# Perturbed angles of attack
alpha_plus = my_alpha + delta_alpha
alpha_minus = my_alpha - delta_alpha

# Run simulations at original and perturbed angles of attack
CD0, CY0, CL0, Cl0, Cm0, Cn0 = get_coefficients(my_alpha)
CD_p, CY_p, CL_p, Cl_p, Cm_p, Cn_p = get_coefficients(alpha_plus)
CD_m, CY_m, CL_m, Cl_m, Cm_m, Cn_m = get_coefficients(alpha_minus)

# Calculate the lift curve slope (C_L_alpha)
C_D_alpha = (CD_p - CD_m) / (2 * delta_alpha)
C_Y_alpha = (CY_p - CY_m) / (2 * delta_alpha)
C_L_alpha = (CL_p - CL_m) / (2 * delta_alpha)
C_l_alpha = (Cl_p - Cl_m) / (2 * delta_alpha)
C_m_alpha = (Cm_p - Cm_m) / (2 * delta_alpha)
C_n_alpha = (Cn_p - Cn_m) / (2 * delta_alpha)

=#

dCF, dCM = stability_derivatives(system)

CDa, CYa, CLa = dCF.alpha
Cla, Cma, Cna = dCM.alpha  # for horizontal tail
CDb, CYb, CLb = dCF.beta   # for vertical tail maybe
Clb, Cmb, Cnb = dCM.beta 
CDp, CYp, CLp = dCF.p
Clp, Cmp, Cnp = dCM.p
CDq, CYq, CLq = dCF.q
Clq, Cmq, Cnq = dCM.q
CDr, CYr, CLr = dCF.r
Clr, Cmr, Cnr = dCM.r

properties = get_surface_properties(system)

Vh= (dh*2*xle_h[2]*yle_h[2])/(2*xle[2]*yle[2]*((chord[1]+chord[2])/2))
Vv= (dv*xle_v[2]*zle_v[2])/(2*xle[2]*yle[2]*2*yle[2])

println("Vh=", Vh)
println("Vv=", Vv)

println("Stability Derivative Pitch: ", "Cma= ", Cma,  " Cmb= ", Cmb,  " Cmp= ", Cmp,  " Cmq= ", Cmq,  " Cmr= ", Cmr)
println("Stability Derivative Yaw: ", "Cna= ", Cna,  " Cnb= ", Cnb,  " Cnp= ", Cnp,  " Cnq= ", Cnq,  " Cnr= ", Cnr)