# Run code for the total number of simulations = 10 000
# Create vectors to store the differences between obs and simulated data to plot
epsilon_mean_r = Float64[]
epsilon_si = Float64[]
epsilon_sinuosity = Float64[]

simulations = 1000
sims = zeros(simulations)
for i = 1:length(sims)

    # Run code for simulated data and "observed" data 10x and get summary stats
    # Create vectors to store the 10x values of summary statistics of sim and obs
    ten_r_s = Float64[]
    ten_si_s = Float64[]
    ten_sinuosity_s = Float64[]
    ten_r_o = Float64[]
    ten_si_o = Float64[]
    ten_sinuosity_o = Float64[]

    iters = 10
    walkers = zeros(iters)
    for i = 1:length(walkers)

        ############################## SIMULATED DATA ##############################

        # Initialize vectors to store the xyz coordinates
        nsteps = 100
        x_s = zeros(nsteps)
        y_s = zeros(nsteps)
        z_s = zeros(nsteps)

        # Set initial time = 0
        t_s = 0

        # Create vectors to store variables
        all_x_s = Float64[]
        all_y_s = Float64[]
        all_z_s = Float64[]
        all_r_s = Float64[]
        turn_angles_s = Float64[]
        time_s = Float64[]

        # Create starting position of the RW at the origin
        x_s[1] = 0.0;
        y_s[1] = 0.0;
        z_s[1] = 0.0;

        # Perform a RW of nsteps
        for i = 2:length(x_s)
            # Sample holding time from exponential distribution or another dist?
            t_next_jump = rand(Exponential())
            # Update the time
            t_s = t_s+t_next_jump

            # Creating a random point in 3D with mean step length = 0.5 and
            # variance = 0.2
            r = rand(TruncatedNormal(0.5,0.1,0,1))
            theta = acos(1-2*rand()) # theta between 0:pi radians
            phi = 2*pi*rand()        # phi between 0:2*pi radians

            # Mapping spherical coordinates onto the cartesian plane
            dx = r*sin(theta)*cos(phi);
            dy = r*sin(theta)*sin(phi);
            dz = r*cos(theta);

            # Updated position
            x_s[i] = x_s[i-1] + dx
            y_s[i] = y_s[i-1] + dy
            z_s[i] = z_s[i-1] + dz

            # Get the current [i] and previous [i-1] coordinates to calculate angle
            # between the 2 vectors = turning angle
            c_1 = x_s[i], y_s[i], z_s[i]
            c_0 = x_s[i-1], y_s[i-1], z_s[i-1]

            # Calculate the turning angle between this vector and previous vector
            turn_angle = acos(vecdot(c_0,c_1)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

            # Push to store all values associated with a coordinate
            push!(all_x_s, x_s[i])
            push!(all_y_s, y_s[i])
            push!(all_z_s, z_s[i])
            push!(all_r_s, r)
            push!(turn_angles_s, turn_angle)
            push!(time_s, t_s)
        end

        ########## SUMMARY STATS SIMULATED DATA #########

        # MEAN STEP LENGTH
        mean_r_s = mean(all_r_s)

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

        sinuosity_s = sd_s / mean_r_s

        # println("mean step length S: ", mean_r_s)
        # println("si S: ", si_s)
        # println("sinuosity S: ", sinuosity_s)

        push!(ten_r_s, mean_r_s)
        push!(ten_si_s, si_s)
        push!(ten_sinuosity_s, sinuosity_s)


        ############################## OBSERVED DATA ##############################

        # Initialize vectors to store the xyz coordinates
        nsteps = 100
        x_o = zeros(nsteps)
        y_o = zeros(nsteps)
        z_o = zeros(nsteps)

        # Set initial time = 0
        t_o = 0

        # Create vectors to store variables
        all_x_o = Float64[]
        all_y_o = Float64[]
        all_z_o = Float64[]
        all_r_o = Float64[]
        turn_angles_o = Float64[]
        time_o = Float64[]

        # Create starting position of the RW at the origin
        x_o[1] = 0.0;
        y_o[1] = 0.0;
        z_o[1] = 0.0;

        # Perform a RW of nsteps
        for i = 2:length(x_o)
            # Sample holding time from exponential distribution or another dist?
            t_next_jump = rand(Exponential())
            # Update the time
            t_o = t_o+t_next_jump

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
            x_o[i] = x_o[i-1] + dx
            y_o[i] = y_o[i-1] + dy
            z_o[i] = z_o[i-1] + dz

            # Get the current [i] and previous [i-1] coordinates to calculate angle
            # between the 2 vectors = turning angle
            c_1 = x_o[i], y_o[i], z_o[i]
            c_0 = x_o[i-1], y_o[i-1], z_o[i-1]

            # Calculate the turning angle between this vector and previous vector
            turn_angle = acos(vecdot(c_0,c_1)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

            # Push to store all values associated with a coordinate
            push!(all_x_o, x_o[i])
            push!(all_y_o, y_o[i])
            push!(all_z_o, z_o[i])
            push!(all_r_o, r)
            push!(turn_angles_o, turn_angle)
            push!(time_o, t_o)
        end

        ########## SUMMARY STATS OBSERVED DATA #########

        # MEAN STEP LENGTH
        mean_r_o = mean(all_r_o)

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

        sinuosity_o = sd_o / mean_r_o

        # println("mean step length o: ", mean_r_o)
        # println("si o: ", si_o)
        # println("sinuosity o: ", sinuosity_o)

        push!(ten_r_o, mean_r_o)
        push!(ten_si_o, si_o)
        push!(ten_sinuosity_o, sinuosity_o)
    end

    # PRINT OUT SS FOR 10X ITERATIONS
    # println("ten_r_s: ", ten_r_s)
    # println("ten_si_s: ", ten_si_s)
    # println("ten_sinuosity_s: ", ten_sinuosity_s)
    # println("**********")
    # println("ten_r_o: ", ten_r_o)
    # println("ten_si_o: ", ten_si_o)
    # println("ten_sinuosity_o: ", ten_sinuosity_o)
    # println("**********")

    # GET AVERAGES FOR THE SS OF 10X ITERATIONS AND PRINT
    average_r_s = mean(ten_r_s)
    average_si_s = mean(ten_si_s)
    average_sinuosity_s = mean(ten_sinuosity_s)

    average_r_o = mean(ten_r_o)
    average_si_o = mean(ten_si_o)
    average_sinuosity_o = mean(ten_sinuosity_o)

    # println("average_r_s: ", average_r_s)
    # println("average_si_s: ", average_si_s)
    # println("average_sinuosity_s: ", average_sinuosity_s)
    # println("**********")
    #
    # println("average_r_o: ", average_r_o)
    # println("average_si_o: ", average_si_o)
    # println("average_sinuosity_o: ", average_sinuosity_o)
    # println("**********")

    # Calculate the difference squared between averages of each summary statistic and
    # push to big epsilon list
    diff_r = (average_r_o - average_r_s)^2
    diff_si = (average_si_o - average_si_s)^2
    diff_sinuosity = (average_sinuosity_o - average_sinuosity_s)^2

    # println("diff_r: ", diff_r)
    # println("diff_si: ", diff_si)
    # println("diff_sinuosity: ", diff_sinuosity)
    # println("**********")

    push!(epsilon_mean_r, diff_r)
    push!(epsilon_si, diff_si)
    push!(epsilon_sinuosity, diff_sinuosity)
end

# PRINT THE EPSILON VECTORS THAT WE WILL PLOT
# Should be the length of simulations
# println("epsilon_mean_r:", epsilon_mean_r)
# println("epsilon_si: ", epsilon_si)
# println("epsilon_sinuosity: ", epsilon_sinuosity)

# PLOT HISTOGRAM OF THE DIFFERENCES OF EACH SUMMARY STATISTIC
# STRAIGHTNESS INDEX
# a = epsilon_si
# plot2 = PyPlot.plt[:hist](a)
# PyPlot.xlabel("Straightness Index Difference between Simulated & Observed")
# PyPlot.title("Straightness Index Difference Distribution")

# SINUOSITY
b = epsilon_sinuosity
plot2 = PyPlot.plt[:hist](b)
PyPlot.xlabel("Epsilon between Simulated & Observed")
PyPlot.title("Sinuosity")
