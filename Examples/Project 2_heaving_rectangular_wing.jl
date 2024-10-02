using VortexLattice
using Plots
pyplot()

# forward velocity
Vinf = 1

# angle of attack
my_alpha = -5*pi/180  # Renamed to avoid conflict

# aspect ratio
AR = 4

# chord lengths
c = [1.0, 0.6, 0.2]

# reduced frequency
k = [0.5, 0.3, 0.1]

t = Vector{Vector{Float64}}(undef, length(k))
CF = Vector{Vector{Vector{Float64}}}(undef, length(k))
CM = Vector{Vector{Vector{Float64}}}(undef, length(k))

for i = 1:length(k)

    # span length
    b = AR*c[i]

    # geometry
    xle = [0.0, 0.0]
    yle = [0.0, b/2]
    zle = [0.0, 0.0]
    chord = [c[i], c[i]]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    ns = 13
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    # reference parameters
    cref = c[i]
    bref = b
    Sref = b*c[i]
    rref = [0.0, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # angular frequency
    ω = 2*Vinf*k[i]/c[i]

    # time
    t[i] = range(0.0, 9*pi/ω, length = 100)
    dt = t[i][2:end] - t[i][1:end-1]
    dt = Vinf*dt

    # heaving amplitude
    h = 0.1*c[i]

    # use forward and vertical velocity at beginning of each time step
    Xdot = Vinf*cos(my_alpha)
    Zdot = Vinf*sin(my_alpha) .- h*cos.(ω*t[i][1:end-1])

    # freestream parameters for each time step
    fs = trajectory_to_freestream(dt; Xdot, Zdot)

    # surface panels
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces = [surface]

    # run analysis
    system, surface_history, property_history, wake_history = unsteady_analysis(
        surfaces, ref, fs, dt; symmetric=symmetric, nwake = 50)

    # extract forces at each time step (uses instantaneous velocity as reference)
    CF[i], CM[i] = body_forces_history(system, surface_history, property_history; frame=Wind())

end

# lift coefficient plot
plot(
    xlim = (6*pi, 8*pi),
    xticks = ([6*pi, 13*pi/2, 7*pi, 15*pi/2, 8*pi], ["\$ 0 \$",
        "\$ \\frac{\\pi}{2} \$", "\$ \\pi \$", "\$ \\frac{3\\pi}{2} \$",
        "\$ 2\\pi \$"]),
    xlabel = "\$ ω \\cdot t \$",
    ylim = (-1.0, 0.1),
    yticks = -1.0:0.2:0.0,
    yflip = true,
    ylabel = "\$ C_{L} \$",
    grid = false,
    )

for i = 1:length(k)
    # extract ω
    ω = 2*Vinf*k[i]/c[i]

    # extract ω*t (use time at the beginning of the time step)
    ωt = ω*t[i][1:end-1]

    # extract CL
    CL = [CF[i][it][3] for it = 1:length(t[i])-1]

    plot!(ωt, CL, label="\$ k = \\frac{\\omega c}{2 U_\\infty} = $(k[i]) \$")
end

plot!(show=true)
