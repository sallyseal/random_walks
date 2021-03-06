# This script does the same as script rw_simulations.jl
# Here we attempt to infer the variance of the theta and phi distributions
# It's an attempt to infer how persistent the random walkl is
# Sample variance from a uniform distribution between 0 and 1.5?
# Ask Michael about range of variance?
# All other parameters such as step length etc. stay constant
# Step length constant at 0.5 mean with variance = 0.1

using Distributions;
using PyPlot;
using StatsBase;
using PyCall, PyPlot; @pyimport seaborn as sns
# Generate the mock data (10x RWs of 100 steps each) and get summary statistics
########## MOCK DATA ##########

# Create vectors to store the average SI and S for 10x RWs
SI_av = Float64[]
S_av = Float64[]

# MSL is kept constant as we are trying to infer variance only
msl = 0.2

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

    # Bounds for distributions
    lower_t = 0
    upper_t = pi
    lower_p = 0
    upper_p = 2*pi

    # Create starting position of the RW at the origin
    x[1] = 0.0;
    y[1] = 0.0;
    z[1] = 0.0;

    # Sample first random point in 3D
    r = rand(TruncatedNormal(msl, 0.1, 0, 1)) # Adam uses log normal?
    theta = acos(1-2*rand()) # theta between 0:pi radians
    phi = 2*pi*rand()        # phi between 0:2*pi radians

    # This is what we want to try to INFER
    # FOR THE PERSISTENCE: variance of theta and phi distributions
    sigma = 0.1 # Can control the tightness/spread of the distribution by altering

    # Perform a RW of nsteps
    for i = 2:length(x)

        # Sample holding time from exponential distribution or another dist?
        t_next_jump = rand(Exponential())
        # Update the time
        t = t+t_next_jump

        # Create variables for updating the distributions
        mu_t = theta
        mu_p = phi

        # Create the distributions for theta and phi to sample next theta and phi
        # Should these be halved?
        # This should be sampled from wrapped normal distribution?
        dist_theta = TruncatedNormal(theta, sigma, lower_t, upper_t)
        dist_phi = TruncatedNormal(phi, sigma, lower_p, upper_p)

        # Randomly sample from the distributions to get updated theta and phi to
        # create next point in 3D space
        theta = rand(dist_theta)
        phi = rand(dist_phi)
        # Here you can change the mean step length we are trying to infer
        r = rand(TruncatedNormal(msl, 0.1, 0, 1))

        # Map spherical point in 3D to the Cartesian Plane
        dx = r*sin(theta)*cos(phi);
        dy = r*sin(theta)*sin(phi);
        dz = r*cos(theta);

        # Updated position
        x[i] = x[i-1] + dx
        y[i] = y[i-1] + dy
        z[i] = z[i-1] + dz

        # Get the coordinate and previous coordinate
        c_0 = x[i], y[i], z[i]
        c_1 = x[i-1], y[i-1], z[i-1]

        # Calculate the angle between this vector and previous vector
        turn_angle = acos(vecdot(c_1,c_0)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

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
variance = Float64[]

# Repeat simulation 10 000x
for i in 1:50000

    # Generate the simulated data (10x RWs of 100 steps each) and get summary stats
    ########## SIMULATED DATA ##########

    # Create vectors to store the average SI and S for 10x RWs
    SI_prime_av = Float64[]
    S_prime_av = Float64[]

    # Sample variance mean from uniform dist between 0 & 2 save value to variance
    v = rand(Uniform(0,2))
    push!(variance, v)

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
        all_theta = Float64[]
        all_phi = Float64[]

        # Bounds for distributions
        lower_t = 0
        upper_t = pi
        lower_p = 0
        upper_p = 2*pi

        # Create starting position of the RW at the origin
        x[1] = 0.0;
        y[1] = 0.0;
        z[1] = 0.0;

        # Sample first random point in 3D
        r = rand(TruncatedNormal(msl, 0.1, 0, 1)) # Adam uses log normal?
        theta = acos(1-2*rand()) # theta between 0:pi radians
        phi = 2*pi*rand()        # phi between 0:2*pi radians

        # FOR THE PERSISTENCE: variance
        sigma = v # Can control the tightness/spread of the distribution by altering

        # Perform a RW of nsteps
        for i = 2:length(x)

            # Sample holding time from exponential distribution or another dist?
            t_next_jump = rand(Exponential())
            # Update the time
            t = t+t_next_jump

            # Create variables for updating the distributions
            mu_t = theta
            mu_p = phi

            # Create the distributions for theta and phi to sample next theta and phi
            # Should these be halved?
            # This should be sampled from wrapped normal distribution?
            dist_theta = TruncatedNormal(theta, sigma, lower_t, upper_t)
            dist_phi = TruncatedNormal(phi, sigma, lower_p, upper_p)

            # Randomly sample from the distributions to get updated theta and phi to
            # create next point in 3D space
            theta = rand(dist_theta)
            phi = rand(dist_phi)
            r = rand(TruncatedNormal(msl, 0.1, 0, 1))

            # Map spherical point in 3D to the Cartesian Plane
            dx = r*sin(theta)*cos(phi);
            dy = r*sin(theta)*sin(phi);
            dz = r*cos(theta);

            # Updated position
            x[i] = x[i-1] + dx
            y[i] = y[i-1] + dy
            z[i] = z[i-1] + dz

            # Get the coordinate and previous coordinate
            c_0 = x[i], y[i], z[i]
            c_1 = x[i-1], y[i-1], z[i-1]

            # Calculate the angle between this vector and previous vector
            turn_angle = acos(vecdot(c_1,c_0)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

            # Push to store all values associated with a coordinate
            push!(all_x, x[i])
            push!(all_y, y[i])
            push!(all_z, z[i])
            push!(all_r, r)
            push!(time, t)
            push!(turn_angles, turn_angle)
            push!(all_theta, theta)
            push!(all_phi, phi)
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

    # Calculate delta and push to delta vector for plotting
    # delta vector will be 10 000 long
    difference_si = sqrt((SI_av - SI_prime_av)^2)
    difference_s = sqrt((S_av - S_prime_av)^2)
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


# CALCULATING THE ACCEPTED V' VALUES FOR PLOTTING
accepted_v_si_1 = Float64[]
accepted_v_si_01 = Float64[]
accepted_v_s_1 = Float64[]
accepted_v_s_01 = Float64[]

zipped_SI = zip(delta_SI, variance)
zipped_S = zip(delta_S, variance)

# PLOTTING THE POSTERIOR DISTRIBUTION OF THE VARIANCE
# Plot the posterior distribution of the variance using S and SI each
# time using 1 and 0.1 percnetiles

# 1. SI_1
for i in zipped_SI
    if i[1] <= e_SI_1
        push!(accepted_v_si_1, i[2])
    end
end
x = accepted_v_si_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: SI_1")

println("size v': ", size(variance))
println("accepted_v_si_1: ", size(accepted_v_si_1))
# ------------------------------------------------------------------------------
# 2. SI_0.1
for i in zipped_SI
    if i[1] <= e_SI_01
        push!(accepted_v_si_01, i[2])
    end
end
x = accepted_v_si_01
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: SI_01")

println("size v': ", size(variance))
println("accepted_v_si_01: ", size(accepted_v_si_01))
# ------------------------------------------------------------------------------
# 3. S_1
for i in zipped_S
    if i[1] <= e_S_1
        push!(accepted_v_s_1, i[2])
    end
end
x = accepted_v_s_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: S_1")

println("size v': ", size(variance))
println("accepted_v_s_1: ", size(accepted_v_s_1))
# ------------------------------------------------------------------------------
# 4. S_0.1
for i in zipped_S
    if i[1] <= e_S_01
        push!(accepted_v_s_01, i[2])
    end
end
x = accepted_v_s_01
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: S_01")

println("size v': ", size(variance))
println("accepted_v_s_01: ", size(accepted_v_s_01))
# ------------------------------------------------------------------------------
