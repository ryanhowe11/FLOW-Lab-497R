using Plots
using VortexLattice
using QuadGK
using Interpolations
using Trapz
#=
function elliptic_wing(root, span, num_sec, filename)
    a = root / 2                # Root length
    b = span                    # Span of a half set of wings
    num_lines = num_sec         # Make formula useful

    # Parametric equations of the half ellipse: x = a * cos(t), y = b * sin(t)
    θ = range(0, π, length=100)
    x = a * cos.(θ)
    y = b * sin.(θ)

    #plot(x, y, aspect_ratio=:equal, legend = false)
    xlabel!("root")
    ylabel!("span")

    # Calculate the y-coordinates of the equidistant lines
    y_lines = range(0, span, length=num_lines+1)

    intersection_points = []
    chord_lengths = Float64[]
    x_intersection_points = Float64[]
    y_intersection_points = Float64[]

    for y_line in y_lines
        x_intersect = a * sqrt(1 - (y_line / b)^2)
        push!(intersection_points, (-x_intersect, y_line))
        push!(intersection_points, (x_intersect, y_line))
        
        #Update X data points
        push!(x_intersection_points, Float64(-x_intersect+a))

        #Update Y data points
        push!(y_intersection_points, Float64(y_line))

        #Update the Chord Lengths Vector
        push!(chord_lengths, Float64(2*x_intersect))

        # Plot the equidistant lines
        #plot!([-x_intersect, x_intersect], [y_line, y_line])
    end

    Area = (pi*a*b)
    cref= 4/(3*pi)*root
    bref=span
    # Save the plot to a PDF file
    #savefig(filename)
    return x_intersection_points, y_intersection_points, intersection_points, chord_lengths, Area, cref, bref
end
=#

"""
    elliptic_wing(root, span, num_sec, filename)
Create and section off an eliptic wing based on its root chord, span, and number of sections
**Arguments**
- `root::Root Chord`: Root Chord of the Wing
- `span::Wing Span`: Span of one wing in isolation
- `num_sec::Number of Sections`: How many sections you want the wing to be discretized into
- `filename::Filename`: File name for a plot of the resulting sectioned wing 
- **Returns**
- `x_intersection_points::X locations of leading edge sections`: The x coordinate where each section starts on the leading edge
- `y_intersection_points::Y locations of leading edge sections`: The y coordinate where each section starts on the leading edge
- `intersection_points::Coordinates of sections on the leading edge`: Where each section starts on the leading edge
- `chord_lengths::Length of section chords`: The chord length at each sectioned divide
- `Area::Reference Area`: The reference area of the wing
"""

L=.5*p*1^2*A*CL

# Define inputs of function
span = 4 #one wing or the whole span       
root = x1
filename = "eliptic_wing_section_plot.pdf"

# geometry (right half of the wing)
xle = (x1, x2, x3,x4, x5, x6, x7)
yle = (y1, y2, y3, y4, y5, y6, y7)
zle = (0, 0, 0, 0, 0, 0 ,0)
chord = (c1, c2, c3, c4, c5, c6, c7)
theta = fill(5*pi/180, num_sec+1)
phi = zeros(num_sec+1)
fc = fill((xc) -> 0, num_sec+1)                     # camberline function for each section
println(xle)
println(yle)
# discretization parameters
ns = 7
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
rref = [0.50, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)
rho = 1.225

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

grids=[grid]

r, c = lifting_line_geometry(grids, 0.25)

cf, cm = lifting_line_coefficients(system, r, c; frame=Wind())

chord=[3, 2.75, 2.5, 2.25, 2, 1.]

for i = 1:ns-1
A+=((chord[i]+chord[i+1])/2)*(yle[i+1]-yle[i])
end

L=.5*rho*Vinf^2*A*CL
