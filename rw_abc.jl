# Script follows same process as rw_simulations, but once we have the list of
# m' values and associated delta_S and delta_SI we iterate over the delta list
# and m' list with zip and if delta <= epsilon then we push the associated m'
# value to the list accepted_m which is the posterior distribution

# We test epsilon for the 1 and 0.1 percentiles

using Distributions;
using PyPlot;
using StatsBase;
using PyCall, PyPlot; @pyimport seaborn as sns
# Generate the mock data (10x RWs of 100 steps each) and get summary statistics
########## MOCK DATA ##########

# Create vectors to store the average SI and S for 10x RWs
SI_av = Float64[]
S_av = Float64[]

random_walks = 10
walks = zeros(random_walks)
for i = 1:length(walks)

    # Initialize vectors to store the xyz coordinates the size of nsteps
    nsteps = 100
    x = zeros(nsteps)
    y = zeros(nsteps)
    z = zeros(nsteps)

    # Set initial time = 0
    t = 0

    # Create vectors to store variables
    all_x = Float64[]
    all_y = Float64[]
    all_z = Float64[]
    all_r = Float64[]
    time = Float64[]
    turn_angles = Float64[]

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

        # Creating a random point in 3D
        r = rand(TruncatedNormal(0.7, 0.1, 0, 1))
        theta = acos(1-2*rand())                # theta between 0:pi radians
        phi = 2*pi*rand()                       # phi between 0:2*pi radians

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
        push!(all_x, x[i])
        push!(all_y, y[i])
        push!(all_z, z[i])
        push!(all_r, r)
        push!(time, t)
        push!(turn_angles, turn_angle)
    end

    # Calculate mock summary statistics

    # Straightness Index: D/L where D= max displacement & L = total path length
    # D = r - r' = sqrt((x-x')^2 + (y-y')^2 + (x-x')^2)
    x1 = all_x[1]
    x2 = all_x[end]
    y1 = all_y[1]
    y2 = all_y[end]
    z1 = all_z[1]
    z2 = all_z[end]
    L = sum(all_r)

    disp = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
    disp = sqrt(disp)
    si = disp/L

    # Sinuosity Index: measures path deviation locally s prop sd/mur
    # where sd = standard dev of turn angle distribution
    # mur = mean step length
    mur = mean(all_r)
    sd = std(turn_angles[2:end])
    s = sd/mur

    # Push ss to vector that stores ss for each one of the 10 runs
    push!(SI_av, si)
    push!(S_av, s)
end
SI_av = mean(SI_av)
S_av = mean(S_av)
println("mock data SI_av: ", SI_av)
println("mock data S_av: ", S_av)



######### SIMULATION 10 000 x ##########

# Create vectors to store deltas for summary stats and mean values used to gen ss
delta_SI = Float64[]
delta_S = Float64[]
means = Float64[]

# Repeat simulation 10 000x
for i in 1:100000

    # Generate the simulated data (10x RWs of 100 steps each) and get summary stats
    ########## SIMULATED DATA ##########

    # Create vectors to store the average SI and S for 10x RWs
    SI_prime_av = Float64[]
    S_prime_av = Float64[]

    # Sample step length mean from uniform dist between 0 & 1 save value to means
    m = rand()
    push!(means, m)

    random_walks = 10
    walks = zeros(random_walks)
    for i = 1:length(walks)

        # Initialize vectors to store the xyz coordinates the size of nsteps
        nsteps = 100
        x = zeros(nsteps)
        y = zeros(nsteps)
        z = zeros(nsteps)

        # Set initial time = 0
        t = 0

        # Create vectors to store variables
        all_x = Float64[]
        all_y = Float64[]
        all_z = Float64[]
        all_r = Float64[]
        time = Float64[]
        turn_angles = Float64[]

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

            # Creating a random point in 3D
            r = rand(TruncatedNormal(m, 0.1, 0, 1))
            theta = acos(1-2*rand())                # theta between 0:pi radians
            phi = 2*pi*rand()                       # phi between 0:2*pi radians

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
            push!(all_x, x[i])
            push!(all_y, y[i])
            push!(all_z, z[i])
            push!(all_r, r)
            push!(time, t)
            push!(turn_angles, turn_angle)
        end

        # Calculate simulated summary statistics

        # Straightness Index: D/L where D= max displacement & L = total path length
        # D = r - r' = sqrt((x-x')^2 + (y-y')^2 + (x-x')^2)
        x1 = all_x[1]
        x2 = all_x[end]
        y1 = all_y[1]
        y2 = all_y[end]
        z1 = all_z[1]
        z2 = all_z[end]
        L = sum(all_r)

        disp = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
        disp = sqrt(disp)
        si = disp/L

        # Sinuosity Index: measures path deviation locally s prop sd/mur
        # where sd = standard dev of turn angle distribution
        # mur = mean step length
        mur = mean(all_r)
        sd = std(turn_angles[2:end])
        s = sd/mur

        # Push ss to vector that stores ss for each one of the 10 runs
        push!(SI_prime_av, si)
        push!(S_prime_av, s)
    end
    SI_prime_av = mean(SI_prime_av)
    S_prime_av = mean(S_prime_av)
    # println(SI_prime_av)
    # println(S_prime_av)

    # Calculate delta and push to delta vector for plotting
    # delta vector will be 10 000 long
    difference_si = (SI_av - SI_prime_av)^2
    difference_s = (S_av - S_prime_av)^2
    # println("difference_si: ", difference_si)
    # println("difference_s: ", difference_s)

    push!(delta_SI, difference_si)
    push!(delta_S, difference_s)
end


# EPSILON CALCULATIONS
# Calculate the 1 and 0.1 percentile of SI and S to generate the epsilon values
e_SI_1 = percentile(delta_SI, 1)
println("e_SI_1: ", e_SI_1)

e_SI_01 = percentile(delta_SI, 0.1)
println("e_SI_01: ", e_SI_01)

e_S_1 = percentile(delta_S, 1)
println("e_S_1: ", e_S_1)

e_S_01 = percentile(delta_S, 0.1)
println("e_S_01: ", e_S_01)


# CALCULATING THE ACCEPTED M' VALUES FOR PLOTTING
accepted_m = Float64[]

zipped_SI = zip(delta_SI, means)
zipped_S = zip(delta_S, means)

# PLOTTING THE POSTERIOR DISTRIBUTION OF THE MEAN STEP LENGTH
# Plot the posterior distribution of the mean step length using S and SI each
# time using 1 and 0.1 percnetiles

# 1. SI_1
# for i in zipped_SI
#     if i[1] <= e_SI_1
#         push!(accepted_m, i[2])
#     end
# end
# x = accepted_m
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length", color="salmon")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution")

# 2. SI_0.1
# for i in zipped_SI
#     if i[1] <= e_SI_01
#         push!(accepted_m, i[2])
#     end
# end
# x = accepted_m
# plot1 = PyPlot.plt[:hist](x; bins=50, alpha=0.5)
# PyPlot.xlabel("Mean Step Length")
# PyPlot.ylabel("Density")
# PyPlot.title("Mean Step Length Posterior Distribution: SI: e = 0.1p")

# 3. S_1
for i in zipped_S
    if i[1] <= e_S_1
        push!(accepted_m, i[2])
    end
end
x = accepted_m
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length", color="salmon")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: S_1")

# 4. S_0.1
# for i in zipped_S
#     if i[1] <= e_S_01
#         push!(accepted_m, i[2])
#     end
# end
# x = accepted_m
# plot1 = PyPlot.plt[:hist](x; bins=50, alpha=0.5)
# PyPlot.xlabel("Mean Step Length")
# PyPlot.ylabel("Density")
# PyPlot.title("Mean Step Length Posterior Distribution: S: e = 0.1p")

println("size m': ", size(means))
println("size accepted_m: ", size(accepted_m))
