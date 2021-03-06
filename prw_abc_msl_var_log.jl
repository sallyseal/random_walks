# This script does the same as script rw_simulations.jl
# Here we attempt to infer mean step length
# We can also attempt to infer how persistent a random walk is - i.e. we can
# attempt to infer the variance of the theta and phi distributions, respectively

using Distributions;
using PyPlot;
using StatsBase;
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd

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

    # Bounds for distributions
    lower_t = 0
    upper_t = pi
    lower_p = 0
    upper_p = 2*pi

    # Create starting position of the RW at the origin
    x[1] = 0.0;
    y[1] = 0.0;
    z[1] = 0.0;

    # PARAMETER 1 TO INFER: msl
    msl = 0.2

    # Sample first random point in 3D
    r = rand(TruncatedNormal(msl, 0.1, 0, 1)) # Adam uses log normal?
    theta = acos(1-2*rand()) # theta between 0:pi radians
    phi = 2*pi*rand()        # phi between 0:2*pi radians

    # PARAMETER 2 TO INFER: variance
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
a = mean(SI_av)
b = std(SI_av)
SI_av = (SI_av - a) / b
SI_mean = mean(SI_av)

c = mean(S_av)
d = std(S_av)
S_av = (S_av - c) / d
S_mean = mean(S_av)

SS = (SI_mean + S_mean)

println("mock data SI_av: ", SI_av)
println("mock data S_av: ", S_av)
println("mock data SS: ", SS)
break
######### SIMULATION 10 000 x ##########

# Create vectors to store deltas for summary stats and mean values used to gen ss
delta_SS = Float64[]
means = Float64[]
variance = Float64[]

# Repeat simulation 10 000x
for i in 1:10000

    # Generate the simulated data (10x RWs of 100 steps each) and get summary stats
    ########## SIMULATED DATA ##########

    # Create vectors to store the average SI and S for 10x RWs
    SI_prime_av = Float64[]
    S_prime_av = Float64[]

    # Sample step length mean from uniform dist between 0 & 1 save value to means
    # Sample variance mean from uniform dist between 0 & 2 save value to variance
    m = rand()
    push!(means, m)
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
        r = rand(TruncatedNormal(m, 0.1, 0, 1)) # Adam uses log normal?
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
            # Here we insert the randomly sampled mean between 0 and 1
            r = rand(TruncatedNormal(m, 0.1, 0, 1))

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
            a = (vecdot(c_1,c_0)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))
            if a <=1 && a >= -1
                turn_angle = acos(a)
            else
                turn_angle = NaN
            end

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
        sd = std(filter(!isnan, turn_angles[2:end]))
        s = sd/mur

        # Push ss to vector that stores ss for each one of the 10 runs
        push!(SI_prime_av, si)
        push!(S_prime_av, s)
    end
    SI_prime_av = (SI_prime_av - (mean(SI_prime_av))) / std(SI_prime_av)
    SI_prime_av = mean(SI_prime_av)
    S_prime_av = (S_prime_av - (mean(S_prime_av))) / std(S_prime_av)
    S_prime_av = mean(S_prime_av)
    SS_prime = (SI_prime_av + S_prime_av)

    # Calculate delta and push to delta vector for plotting
    # delta vector will be 10 000 long
    difference_ss = sqrt((SS - SS_prime)^2)

    # println("difference_si: ", difference_si)
    # println("difference_s: ", difference_s)
    push!(delta_SS, difference_ss)
end


x = delta_SS
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Delta Combined Summary Statistic")
ax[:set_xlim]
ax[:set_title]("Distribution of Delta Combined Summary Statistic: PRW: SS_1")

# EPSILON CALCULATIONS
# Calculate the 1 and 0.1 percentile of SS to generate the epsilon values
e_SS_1 = percentile(delta_SS, 1)
println("e_SS_1: ", e_SS_1)

e_SS_01 = percentile(delta_SS, 0.1)
println("e_SS_01: ", e_SS_01)

# CALCULATING THE ACCEPTED M' & V' VALUES FOR PLOTTING
accepted_m_1 = Float64[]
accepted_m_01 = Float64[]

accepted_v_1 = Float64[]
accepted_v_01 = Float64[]

zipped_SS = zip(delta_SS, means, variance)

# PLOTTING THE POSTERIOR DISTRIBUTION OF THE MEAN STEP LENGTH
# Plot the posterior distribution of the mean step length using S and SI each
# time using 1 and 0.1 percnetiles
# -----------------------------------------------------------------------------
# 1. SS_1
for i in zipped_SS
    if i[1] <= e_SS_1
        push!(accepted_m_1, i[2])
        push!(accepted_v_1, i[3])
    end
end
x = accepted_m_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution Combined: PRW: SS_1")

y = accepted_v_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution Combined: PRW: SS_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
sns.jointplot(x="msl", y="variance", data=df, kind="kde")

println("size m': ", size(means))
println("size v': ", size(variance))
println("size accepted_m_1: ", size(accepted_m_1))
println("size accepted_v_1: ", size(accepted_v_1))
#------------------------------------------------------------------------------
#
# # 2. SS_01
# for i in zipped_SS
#     if i[1] <= e_SS_01
#         push!(accepted_m_01, i[2])
#         push!(accepted_v_01, i[3])
#     end
# end
# x = accepted_m_01
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution Combined: PRW: SS_0.1")
#
# y = accepted_v_01
# fig,ax = PyPlot.subplots()
# sns.distplot(y, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution Combined: PRW: SS_01")
#
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")
#
# println("size m': ", size(means))
# println("size v': ", size(variance))
# println("size accepted_m_01: ", size(accepted_m_01))
# println("size accepted_v_01: ", size(accepted_v_01))
# ------------------------------------------------------------------------------

# a = sort(accepted_m_01)
# fig,ax = PyPlot.subplots()
# accepted_m_100 = a[1:101]
# sns.distplot(a, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution Combined: PRW: SS_1: Best 100")
#
# b = sort(accepted_v_01)
# fig,ax = PyPlot.subplots()
# accepted_v_100 = b[1:101]
# sns.distplot(b, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution Combined: PRW: SS_1: Best 100")
#
# a = x
# b = y
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")






# SI_prime_av = (SI_prime_av - (mean(SI_prime_av))) / std(SI_prime_av)
# SI_prime_av = log(mean(SI_prime_av) + 1)
# S_prime_av = (S_prime_av - (mean(S_prime_av))) / std(S_prime_av)
# S_prime_av = log(mean(S_prime_av) + 1)
# SS_prime = (SI_prime_av + S_prime_av)
