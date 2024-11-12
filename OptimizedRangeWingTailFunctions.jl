using LinearAlgebra
    
#Setting up all the problem variables
function ProblemSetup(num_sec, density, c)
T = eltype(c)
    # Define inputs of function
 span = 8.0 #one wing or the whole span       
 rho = 1.225
 Vinf = 2*c[num_sec+2]

chords = zeros(num_sec+1)

for i in 1:num_sec+1
chords[i] = c[i]
end

# Calculate the average
average_chord = mean(chords)

dt=c[num_sec+3]
lh=c[num_sec+4]
lv=c[num_sec+5]

weight=1.7*density
# weight = ((average_chord*span+lh*lh/4+lv*lv/4)*.333333+dt*pi*.225)*density


 yle=zeros(T, num_sec+1)
 # geometry (right half of the wing)
 yle = [i * (span / (num_sec)) for i in 0:(num_sec)]
 zle = zeros(T, num_sec+1)
 theta = zeros(T, num_sec+1)
 phi = zeros(T, num_sec+1)
 fc = zeros(T, num_sec+1)
 xle = zeros(T, num_sec+1)

 for i in 1:num_sec
    xle[i+1] = (chords[1]/4 - chords[i+1]/4)
end

     # discretization parameters
 ns = num_sec
 nc = num_sec
 spacing_s = Uniform()
 spacing_c = Uniform()

return lh, lv, dt, chords, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c
end

#Find the leading edge values assuming zero sweep
function XLEcalc(xle, chords, num_sec)

 for i in 1:num_sec
     xle[i+1] = (chords[1]/4 - chords[i+1]/4)
 end

 return xle
end

#Define your references such as area and position
function ReferenceCalculation(num_sec, yle, span, Vinf, c)

 Sref= 0.0

 # Reference Area Calculation
 for i in 1:num_sec
     S = ((c[i] + c[i+1]) / 2) * (yle[i+1] - yle[i])
     Sref=S+Sref
 end

 cref= Sref/span

 # reference parameters
 rref = [0.50, 0.0, 0.0]
 ref = Reference(Sref, cref, span, rref, Vinf)

return Sref, ref
end

#gives you your fs variable to use in Vortex Lattice
function FreestreamParams(Vinf)

 # freestream parameters
 alpha_angle = 5.0*pi/180
 beta = 0.0
 Omega = [0.0; 0.0; 0.0]
 fs = Freestream(Vinf, alpha_angle, beta, Omega)

return fs
end

#Finds the coefficients of lift and drag for whole system
function CalculateSurfaceForcesCoeff(surface, htail, vtail, ref, fs, symmetric)

 # create vector containing all surfaces
 surfaces = [surface, htail, vtail]
 surface_id = [1, 2, 3]

 # perform steady state analysis
 system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

 # retrieve near-field forces
 CF, CM = body_forces(system; frame=Wind())

 CD, CY, CL = CF
 Cl, Cm, Cn = CM

return CL, CD
end

#Setting Up Optimization Problem Variables
function OptimizationSetup(num_sec, xstart)
 # Initialize vectors based on num_sec
 c0 = xstart  # starting point

 # if num_sec >= length(chord_opt)
 # for i in 1:length(chord_opt)
 # c0[i]=chord_opt[i]
 # end
 # else
 #     for i in 1:num_sec
 #     c0[i]=chord_opt[i]
 #     end
 # end

 lc = fill(0.01, num_sec+5)  # lower bounds on x
 uc = fill(5.0, num_sec+5)  # upper bounds on x
 ng = 4 + 2*num_sec  # number of constraints
 lg = -Inf*ones(ng)  # lower bounds on g
 ug = zeros(ng)  # upper bounds on g
 g = zeros(ng)

 # ----- set some options ------
 ip_options = Dict(
     "max_iter" => 500,
     "tol" => 1e-3
 )
 solver = IPOPT(ip_options)
 options = Options(;solver, derivatives=ForwardFD())
 return c0, ng, lc, uc, lg, ug, options, g
end

# Plotting function
function plot_chords(xle_opt, yle, chords, num_sec)

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

function RunOptimizer(num_sec, density, xstart)

c0, ng, lc, uc, lg, ug, options, g = OptimizationSetup(num_sec, xstart)

#Creating the optimization problem
function wing_optimizer(g, c)

    lh, lv, dt, chords, span, rho, Vinf, weight, yle, zle, theta, phi, fc, xle, ns, nc, spacing_s, spacing_c = ProblemSetup(num_sec, density, c)

    xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v=SetUpTail(c, lh, lv)

    # xle = XLEcalc(xle, chords, num_sec)

    Sref, ref = ReferenceCalculation(num_sec, yle, span, Vinf, c)

    fs = FreestreamParams(Vinf)

    # construct surface
    grid, surface = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
    spacing_s=spacing_s, spacing_c=spacing_c)

    # generate surface panels for horizontal tail
    hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    VortexLattice.translate!(hgrid, [dt, 0.0, 0.0])
    VortexLattice.translate!(htail, [dt, 0.0, 0.0])

    # generate surface panels for vertical tail
    vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    VortexLattice.translate!(vgrid, [dt, 0.0, 0.0])
    VortexLattice.translate!(vtail, [dt, 0.0, 0.0])

        # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = [true, true, false]

    CL, CD = CalculateSurfaceForcesCoeff(surface, htail, vtail, ref, fs, symmetric)

    Cla, Cma, Cmq, CYb, CLb, Cnb, CLp, Cnr = GetStabilityDerivatives(surface, htail, vtail, ref, fs, symmetric)

    D=.5*rho*Vinf^2*Sref*CD
    L = .5*rho*Vinf^2*Sref*CL

    g[1]= weight - L
    g[2]= c[1] - c[2] - 0.02
    g[3]= Cma-0.1
    # g[4]= Cmq-0.1
    g[4]= -CYb
    # g[6]= Cnr-0.1
    # g[3]= -Cla
    # g[4]= Cma
    # g[5]= Cmq
    # g[6]= -CYb
    # g[7]= CLb
    # g[8]= -Cnb
    # g[9]= CLp
    # g[10]= Cnr

# CLa +, CMa -, CMq -, CYb -, CLb -, CNb +, CLp -, CNr -

    # Calculate chord differences
    for i in 1:num_sec
        g[i+4] = c[i+1] - c[i]
    end

        # Calculate chord differences
    for i in 1:num_sec
        g[i+num_sec+4] = c[i] - c[i+1]-.5
    end

    O=D/(sqrt(abs(L))*Vinf)

return O
end

xopt, fopt, info = minimize(wing_optimizer, c0, ng, lc, uc, lg, ug, options)

return xopt, fopt, info
end

function SetUpTail(c, lh, lv)

    T = eltype(c)
    # horizontal stabilizer
    xle_h = Array{T}([0.0, lh/4])
    yle_h = Array{T}([0.0, lh])
    zle_h = Array{T}([0.1, 0.1])
    chord_h = Array{T}([0.7, 0.42])
    theta_h = Array{T}([0.0, 0.0])
    phi_h = Array{T}([0.0, 0.0])
    fc_h = fill((xc) -> 0, 2) #camberline function for each section
    ns_h = 6
    nc_h = 3
    spacing_s_h = Uniform()
    spacing_c_h = Uniform()
    mirror_h = false

    # vertical stabilizer
    xle_v = Array{T}([0.0, lv/4])
    yle_v = Array{T}([0.0, 0.0])
    zle_v = Array{T}([0.1, lv])
    chord_v = Array{T}([0.7, 0.42])
    theta_v = Array{T}([0.0, 0.0])
    phi_v = Array{T}([0.0, 0.0])
    fc_v = fill((xc) -> 0, 2) #camberline function for each section
    ns_v = 5
    nc_v = 3
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

return xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h, spacing_s_h, spacing_c_h, mirror_h, xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v, spacing_s_v, spacing_c_v, mirror_v
end

function GetStabilityDerivatives(surface, htail, vtail, ref, fs, symmetric)
# create vector containing all surfaces
surfaces = [surface, htail, vtail]
surface_id = [1, 2, 3]

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)


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
return Cla, Cma, Cmq, CYb, CLb, Cnb, CLp, Cnr
end

# Put a constraint 