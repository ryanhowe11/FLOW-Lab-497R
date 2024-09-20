# Katz and Plotkin: Figures 13.34 and 13.35
# AR = [4, 8, 12, 20, ∞]
# Vinf*Δt/c = 1/16
# α = 5°

using VortexLattice
using Plots
pyplot()

AR = [4, 8, 12, 20, 1e3] # last aspect ratio is essentially infinite

system = Vector{Any}(undef, length(AR))
surface_history = Vector{Any}(undef, length(AR))
property_history = Vector{Any}(undef, length(AR))
wake_history = Vector{Any}(undef, length(AR))
CF = Vector{Vector{Vector{Float64}}}(undef, length(AR))
CM = Vector{Vector{Vector{Float64}}}(undef, length(AR))

# non-dimensional time (t*Vinf/c)
t = range(0.0, 10.0, step=1/16)

# chord length
c = 1

# time step
dt = [t[i+1]-t[i] for i = 1:length(t)-1]

for i = 1:length(AR)

    # span length
    b = AR[i]*c

    # planform area
    S = b*c

    # geometry
    xle = [0.0, 0.0]
    yle = [-b/2, b/2]
    zle = [0.0, 0.0]
    chord = [c, c]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    ns = 13
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = false

    # reference parameters
    cref = c
    bref = b
    Sref = S
    rref = [0.0, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # create vortex rings
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing surfaces
    surfaces = [surface]

    # run analysis
    system[i], surface_history[i], property_history[i], wake_history[i] =
        unsteady_analysis(surfaces, ref, fs, dt; symmetric, wake_finite_core = false)

    # extract forces at each time step
    CF[i], CM[i] = body_forces_history(system[i], surface_history[i],
        property_history[i]; frame=Wind())

end

write_vtk("acceleration-AR4", surface_history[1], property_history[1],
    wake_history[1], dt; symmetric=false)


# lift coefficient plot
plot(
    xlim = (0.0, 10.0),
    xticks = 0.0:1.0:10.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 0.55),
    yticks = 0.0:0.1:0.5,
    ylabel = "\$ C_{L} \$",
    grid = false,
    overwrite_figure=false
    )

for i = 1:length(AR)
    CL = [CF[i][j][3] for j = 1:length(CF[i])]
    plot!(t[2:end], CL, label="AR = $(AR[i])")
end

plot!(show=true)

# drag coefficient plot
plot(
    xlim = (0.0, 10.0),
    xticks = 0.0:1.0:10.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 0.030),
    yticks = 0.0:0.005:0.03,
    ylabel = "\$ C_{D} \$",
    grid = false,
    overwrite_figure=false
    )

for i = 1:length(AR)
    CD = [CF[i][j][1] for j = 1:length(CF[i])]
    plot!(t[2:end], CD, label="AR = $(AR[i])")
end

plot!(show=true)