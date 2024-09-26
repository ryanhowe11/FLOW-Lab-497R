using Plots
using VortexLattice

function elliptic_wing(root, span, num_sec, filename)
    a = root / 2                # Root length
    b = span                    # Span of a half set of wings
    num_lines = num_sec         # Make formula useful

    # Parametric equations of the half ellipse: x = a * cos(t), y = b * sin(t)
    θ = range(0, π, length=100)
    x = a * cos.(θ)
    y = b * sin.(θ)

    plot(x, y, aspect_ratio=:equal, legend = false)
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
        push!(x_intersection_points, Float64(x_intersect))

        #Update Y data points
        push!(y_intersection_points, Float64(y_line))

        #Update the Chord Lengths Vector
        push!(chord_lengths, Float64(2*x_intersect))

        # Plot the equidistant lines
        plot!([-x_intersect, x_intersect], [y_line, y_line])
    end



    # Save the plot to a PDF file
    savefig(filename)
    return x_intersection_points, y_intersection_points, intersection_points, chord_lengths
end

for i in 1:40
# Define inputs of function
span = 10       
root = 5
num_sec = i
filename = "eliptic_wing_section_plot.pdf"
x_points, y_points, points, chords = elliptic_wing(root, span, num_sec, filename)
#println("Intersection points: ", points)
#println("Section Chord Lengths: ", chords)
#println("Plot saved to ", filename)

# geometry (right half of the wing)
xle = x_points
yle = y_points
zle = zeros(num_sec+1)
chord = chords
theta = fill(2.0*pi/180, num_sec+1)
phi = zeros(num_sec+1)
fc = fill((xc) -> 0, num_sec+1)                     # camberline function for each section

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
alpha_angle = 1.0*pi/180
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

properties = get_surface_properties(system)

#write_vtk("symmetric-planar-wing", surfaces, properties; symmetric)

# construct geometry with mirror image
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc=fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

# symmetry is not used in the analysis
symmetric = false

# create vector containing all surfaces
surfaces = [surface]

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

# Calculate aerodynamic efficiency
efficiency = CL / CD
println("Aerodynamic Efficiency (L/D ratio): ", efficiency, " Number of Sections: ", num_sec)
#write_vtk("mirrored-planar-wing", surfaces, properties; symmetric)
end