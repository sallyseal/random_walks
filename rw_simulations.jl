# This script...

# Import necessary libraries
using Gadfly;
using Distributions;
using PyPlot;
import Plots;

# Create vector to store differences between summary statistics for sim vs. "observed"
# Will plot the distribution of these differences at end
epsilon_mean_r = zeros(iterations)
epsilon_si = zeros(iterations)
epsilon_sinuosity = zeros(iterations)

# Number of iterations to complete = 10 000
iterations = 5

for i = 1:length(iterations)

    # Create vectors to store the summary stats for simulated and "observed"
    # Will calculate the differences in simulated vs observed and push these
    # differences to epsilon
    sim_r_mean = Float64[]
    sim_si = Float64[]
    sim_sinuosity = Float64[]

    observed_r_mean = Float64[]
    observed_si = Float64[]
    observed_sinuosity = Float64[]

    # Run code for simulated data and "observed" data 10x and get summary stats
    iters = 10

    for i = 1:length(iters)

###############################################################################

        ### SIMULATED DATA ###

        # Initialize vectors to store the xyz coordinates
        nsteps = 100
        x_s = zeros(nsteps)
        y_s = zeros(nsteps)
        z_s = zeros(nsteps)

        # Set initial time = 0
        t = 0

        # Create vectors to store variables
        all_x_s = Float64[]
        all_y_s = Float64[]
        all_z_s = Float64[]
        all_r_s = Float64[]
        turn_angles_s = Float64[]
        time_s = Float64[]

        # Create starting position of the RW at the origin
        x[1] = 0.0;
        y[1] = 0.0;
        z[1] = 0.0;

        # Perform a RW of nsteps
        for i = 2:length(x)
            # Sample holding time from exponential distribution or another dist?
            t_next_jump = rand(Exponential())
            # Update the time
            t = t+t_next_jump

            # Creating a random point in 3D with mean step length = 0.5 and
            # variance = 0.2
            r = rand(TruncatedNormal(0.5,0.2,0,1))
            theta = acos(1-2*rand()) # theta between 0:pi radians
            phi = 2*pi*rand()        # phi between 0:2*pi radians

            # Mapping spherical coordinates onto the cartesian plane
            dx = r*sin(theta)*cos(phi);
            dy = r*sin(theta)*sin(phi);
            dz = r*cos(theta);

            # Updated position
            x[i] = x[i-1] + dx
            y[i] = y[i-1] + dy
            z[i] = z[i-1] + dz

            # Get the current [i] and previous [i-1] coordinates to calculate angle
            # between the 2 vectors = turning angle
            c_1 = x[i], y[i], z[i]
            c_0 = x[i-1], y[i-1], z[i-1]

            # Calculate the turning angle between this vector and previous vector
            turn_angle = acos(vecdot(c_0,c_1)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

            # Push to store all values associated with a coordinate
            push!(all_x_s, x[i])
            push!(all_y_s, y[i])
            push!(all_z_s, z[i])
            push!(all_r_s, r)
            push!(turn_angles_s, turn_angle)
            push!(time_s, t)
        end

        # Calculate summary statistics for the 10x simulated runs and take average
        # Create vectors to store averages
        mean_r_s = Float64[]
        mean_si_s = Float64[]
        mean_sinuosity_s = Float64[]

        # MEAN STEP LENGTH
        mean_r = mean(all_r_s)

        # STRAIGHTNESS INDEX: D/L
        x1 = all_x_s[1]
        x2 = all_x_s[end]
        y1 = all_y_s[1]
        y2 = all_y_s[end]
        z1 = all_z_s[1]
        z2 = all_z_s[end]
        D_s = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
        D_s = sqrt(D_s)
        L = sum(all_r_s)

        si_s = D_s / L

        # SINUOSITY: sd of turn angles / mean step length
        sd_s = std(turn_angles_s[2:end])
        sinuosity_s = sd_s / mean_r

        # Push the average of the summary statistics for 10x runs to a vector
        push!(mean_r_s, mean_r)
        push!(mean_si_s, mean(si_s))
        push!(mean_sinuosity_s, mean(sinuosity_s))

###############################################################################

        ### OBSERVED DATA ###

        # Initialize vectors to store the xyz coordinates
        nsteps = 100
        x_o = zeros(nsteps)
        y_o = zeros(nsteps)
        z_o = zeros(nsteps)

        # Set initial time = 0
        t = 0

        # Create vectors to store variables
        all_x_o = Float64[]
        all_y_o = Float64[]
        all_z_o = Float64[]
        all_r_o = Float64[]
        turn_angles_o = Float64[]
        time_o = Float64[]

        # Create starting position of the RW at the origin
        x[1] = 0.0;
        y[1] = 0.0;
        z[1] = 0.0;

        # Perform a RW of nsteps
        for i = 2:length(x)
            # Sample holding time from exponential distribution or another dist?
            t_next_jump = rand(Exponential())
            # Update the time
            t = t+t_next_jump

            # Creating a random point in 3D with mean step length = 0.8 and
            # variance = 0.2
            r = rand(TruncatedNormal(0.8,0.2,0,1))
            theta = acos(1-2*rand()) # theta between 0:pi radians
            phi = 2*pi*rand()        # phi between 0:2*pi radians

            # Mapping spherical coordinates onto the cartesian plane
            dx = r*sin(theta)*cos(phi);
            dy = r*sin(theta)*sin(phi);
            dz = r*cos(theta);

            # Updated position
            x[i] = x[i-1] + dx
            y[i] = y[i-1] + dy
            z[i] = z[i-1] + dz

            # Get the current [i] and previous [i-1] coordinates to calculate angle
            # between the 2 vectors = turning angle
            c_1 = x[i], y[i], z[i]
            c_0 = x[i-1], y[i-1], z[i-1]

            # Calculate the turning angle between this vector and previous vector
            turn_angle = acos(vecdot(c_0,c_1)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

            # Push to store all values associated with a coordinate
            push!(all_x_o, x[i])
            push!(all_y_o, y[i])
            push!(all_z_o, z[i])
            push!(all_r_o, r)
            push!(turn_angles_o, turn_angle)
            push!(time_o, t)
        end

        # Calculate summary statistics for the 10x simulated runs and take average
        # Create vectors to store averages
        mean_r_o = Float64[]
        mean_si_o = Float64[]
        mean_sinuosity_o = Float64[]

        # MEAN STEP LENGTH
        mean_r = mean(all_r_o)

        # STRAIGHTNESS INDEX: D/L
        x1 = all_x_o[1]
        x2 = all_x_o[end]
        y1 = all_y_o[1]
        y2 = all_y_o[end]
        z1 = all_z_o[1]
        z2 = all_z_o[end]
        D_o = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
        D_o = sqrt(D_o)
        L = sum(all_r_o)

        si_o = D_o / L

        # SINUOSITY: sd of turn angles / mean step length
        sd_o = std(turn_angles_o[2:end])
        sinuosity_o = sd_o / mean_r

        # Push the average of the summary statistics for 10x runs to a vector
        push!(mean_r_o, mean_r)
        push!(mean_si_o, mean(si_o))
        push!(mean_sinuosity_o, mean(sinuosity_o))

###############################################################################
    end

    # Get the difference between the average of ss between simulated and observed
    # for the 10x iterations. Do this 10 000 times and push these 10 000 differences
    # to epsilon vectors for plotting
    # Can take square of difference or absolute value of the difference?

    diff_r = (mean_r_o - mean_r_s)^2
    diff_si = (mean_si_o - mean_si_s)^2
    diff_sinuosity = (mean_sinuosity_o - mean_sinuosity_s)^2

    push!(epsilon_mean_r, diff_r)
    push!(epsilon_si, diff_si)
    push!(epsilon_sinuosity, diff_sinuosity)

# CALCULATE THE DIFFERENCES BETWEEN SIMULATED AND OBSERVED AND PUSH TO
# EPSILON VECTORS FOR PLOTTING
println("epsilon_mean_r: ", epsilon_mean_r)
println("epsilon_si: ", epsilon_si)
println("epsilon_sinuosity: ", epsilon_sinuosity)
