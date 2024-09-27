using Plots
using VortexLattice

function tapered_wing(root, span, tip, num_sec, filename)
    a = root / 2                # Root length
    b = span                    # Span of a half set of wings
    c = tip / 2
    num_lines = num_sec         # Make formula useful

    # Coordinates of the trapezoid
    x = [-a, a, c, -c, -a]
    y = [0, 0, b, b, 0]

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
        # Linear interpolation between root and tip
        x_intersect = a - (a - c) * (y_line / b)
        push!(intersection_points, (-x_intersect, y_line))
        push!(intersection_points, (x_intersect, y_line))
        
        #Update X data points
        push!(x_intersection_points, Float64(-x_intersect+a))

        #Update Y data points
        push!(y_intersection_points, Float64(y_line))

        #Update the Chord Lengths Vector
        push!(chord_lengths, Float64(2*x_intersect))

        # Plot the equidistant lines
        plot!([-x_intersect, x_intersect], [y_line, y_line])
    end



    # Save the plot to a PDF file
    savefig(filename)

    Area = span*(tip+root)/2
    cref= (tip+root)/2
    bref=2*span
    # Save the plot to a PDF file
    savefig(filename)
    return x_intersection_points, y_intersection_points, intersection_points, chord_lengths, Area, cref, bref
end

"""
    elliptic_wing(root, span, num_sec, filename)
Create and section off an eliptic wing based on its root chord, span, and number of sections
**Arguments**
- `root::Root Chord`: Root Chord of the Wing
- `span::Wing Span`: Span of one wing in isolation
- `tip::Tip Chord`: Tip Chord of the Wing
- `num_sec::Number of Sections`: How many sections you want the wing to be discretized into
- `filename::Filename`: File name for a plot of the resulting sectioned wing 
- **Returns**
- `x_intersection_points::X locations of leading edge sections`: The x coordinate where each section starts on the leading edge
- `y_intersection_points::Y locations of leading edge sections`: The y coordinate where each section starts on the leading edge
- `intersection_points::Coordinates of sections on the leading edge`: Where each section starts on the leading edge
- `chord_lengths::Length of section chords`: The chord length at each sectioned divide
- `Area::Reference Area`: The reference area of the wing
"""

# Define inputs of function
span = 10       
root = 3
tip = 1.5
num_sec = 20
filename = "tapered_wing_section_plot.pdf"
x_points, y_points, points, chords, Sref, cref, bref = tapered_wing(root, span, tip, num_sec, filename)
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
ns = length(yle)
nc = 1
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
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

# Calculate aerodynamic efficiency
#println("Aerodynamic Efficiency (L/D ratio): ", efficiency, " Number of Sections: ", num_sec)
write_vtk("tapered-mirrored-planar-wing", surfaces, properties; symmetric)

# Assuming cf is your vector of matrices
z_direction_coefficients = []

# Loop through each matrix in the cf vector
for matrix in cf
    # Extract the third row (z-direction force coefficients)
    z_coefficients = matrix[3, :]
    push!(z_direction_coefficients, z_coefficients)
end

y = collect(range(0, stop=span, step=0.1))

# Convert the list of z-direction coefficients to an array if needed
lifting_coefficients = hcat(z_direction_coefficients...)

Lift_prime = .5*rho*Vinf^2*lifting_coefficients[:, 1] .* chords

area_prime = trapz(yle, Lift_prime)

bprime = span
aprime = 4*area_prime/(bprime*pi)

# Calculate the ideal elliptic lift distribution
θ = range(0, π/2, length=100)
x = bprime * cos.(θ)
y = aprime * sin.(θ)

# Calculate the ideal elliptic lift distribution
#Cl_max = maximum(Lift_prime)
#elliptical_distribution = Cl_max * sqrt.(1 .- (y ./ span).^2)

# Plot the lift distribution
plot(yle, Lift_prime, label="Calculated Lift Distribution", xlabel="Spanwise Location (y)", ylabel="Lift Coefficient (Cl)")
#plot!(y, elliptical_distribution, label="Elliptical Lift Distribution", linestyle=:dash)
plot!(x, y, linestyle=:dash)

# Save the lift distribution plot as a PDF
savefig("Lift_Distribution_along_the_Span_taper.pdf")