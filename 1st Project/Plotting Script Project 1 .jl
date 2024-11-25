using Xfoil, Plots, Printf, DelimitedFiles, CSV, DataFrames
pyplot()

function read_airfoil_coordinates(filename::String)
    # Read airfoil coordinates from a file
    data = readdlm(filename)
    
    # Separate the columns into x and y
    x = data[:, 1]
    y = data[:, 2]

    # Load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x, y)
    return x, y
end

function repanel_airfoil()
    # Repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane()
    return xr, yr
end

function set_operating_conditions()
    # Set operating conditions
    aoa_range = -9:1:14  # Range of angle of attacks, in degrees
    re = 1.0e5           # Reynolds number
    mach = 0.0           # Mach #
    complex = 0
    gr()
    return aoa_range, re, mach, complex
end

function solve_alpha(aoa_range::StepRange{Int64, Int64}, re::Float64)
    # Initialize outputs
    n_a = length(aoa_range)
    c_l = zeros(n_a)
    c_d = zeros(n_a)
    c_dp = zeros(n_a)
    c_m = zeros(n_a)
    converged = zeros(Bool, n_a)

    # Determine airfoil coefficients across a range of angle of attacks
    for i = 1:n_a
        c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(aoa_range[i], re; iter=100, reinit=true)
    end

    # Print results
    println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
    for i = 1:n_a
        @printf("%8f\t%8f\t%8f\t%8f\t%d\n", aoa_range[i], c_l[i], c_d[i], c_m[i], converged[i])
    end
    return c_l, c_d, c_m
end

function main(filename::String)
    # Read airfoil coordinates
    x, y = read_airfoil_coordinates(filename)
    
    # Set operating conditions
    aoa_range, re, mach, complex = set_operating_conditions()
    
    # Repanel airfoil
    xr, yr = repanel_airfoil()

    # Solve airfoil
    c_l, c_d, c_m = solve_alpha(aoa_range, re)

    return xr, yr, aoa_range, c_l, c_d, c_m, complex
end

# Call the main function and store the results
filename = "naca2412.txt"
airfoil = "Moment New"
xr, yr, aoa_range, c_l, c_d, c_m, complex = main(filename)

# Plotting code

# Read the Third Party Data CSV file
df = CSV.File("0012 Drag.csv", header=false) |> DataFrame

# Plot Third Party Data the data
#plot!(p, df.Column1, df.Column2, label = "Third Party NACA 0012")

# Plot the coefficients
#p = plot(aoa_range, c_l)



p=plot(aoa_range, c_m, grid=false, label =  "Generated NACA 0012", xlabel="Angle of Attack (degrees)", ylabel="Moment Coefficient", show=true, legend=false)
savefig(p, "$(airfoil).pdf")