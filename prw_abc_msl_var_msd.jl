# This script does the same as script rw_simulations.jl
# Here we attempt to infer mean step length
# We can also attempt to infer how persistent a random walk is - i.e. we can
# attempt to infer the variance of the theta and phi distributions, respectively
# combined = SI + MSD

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
msd_av = Float64[]

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

        # Calculate msd
        msd = (c_0[1] - c_1[1])^2 + (c_0[2] - c_1[2])^2 + (c_0[3] - c_1[3])^2
        # msd = sqrt(msd)

        # Push to store all values associated with a coordinate
        push!(all_x, x[i])
        push!(all_y, y[i])
        push!(all_z, z[i])
        push!(all_r, r)
        push!(time, t)
        push!(turn_angles, turn_angle)
        push!(msd_av, msd)
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
msd_av = mean(msd_av)
combined_av = (SI_av + msd_av)
println("mock data SI_av: ", SI_av)
println("mock data S_av: ", S_av)
println("mock data msd_av: ", msd_av)
println("mock data combined: ", combined_av)


######### SIMULATION 10 000 x ##########

# Create vectors to store deltas for summary stats and mean values used to gen ss
delta_SI = Float64[]
delta_S = Float64[]
delta_msd = Float64[]
delta_combined = Float64[]
means = Float64[]
variance = Float64[]

# Repeat simulation 10 000x
for i in 1:1000

    # Generate the simulated data (10x RWs of 100 steps each) and get summary stats
    ########## SIMULATED DATA ##########

    # Create vectors to store the average SI and S for 10x RWs
    SI_prime_av = Float64[]
    S_prime_av = Float64[]
    msd_prime_av = Float64[]

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

            # Calculate msd
            msd = (c_0[1] - c_1[1])^2 + (c_0[2] - c_1[2])^2 + (c_0[3] - c_1[3])^2
            # msd = sqrt(msd)

            # Push to store all values associated with a coordinate
            push!(all_x, x[i])
            push!(all_y, y[i])
            push!(all_z, z[i])
            push!(all_r, r)
            push!(time, t)
            push!(turn_angles, turn_angle)
            push!(msd_prime_av, msd)
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
    SI_prime_av = mean(SI_prime_av)
    S_prime_av = mean(S_prime_av)
    msd_prime_av = mean(msd_prime_av)
    combined_prime_av = mean(SI_prime_av + msd_prime_av)

    # Calculate delta and push to delta vector for plotting
    # delta vector will be 10 000 long
    difference_si = sqrt((SI_av - SI_prime_av)^2)
    difference_s = sqrt((S_av - S_prime_av)^2)
    difference_msd = sqrt((msd_av - msd_prime_av)^2)
    difference_combined = sqrt((combined_av - combined_prime_av)^2)
    # println("difference_si: ", difference_si)
    # println("difference_s: ", difference_s)

    push!(delta_SI, difference_si)
    push!(delta_S, difference_s)
    push!(delta_msd, difference_msd)
    push!(delta_combined, difference_combined)
end
println("delta_SI: ", delta_SI)
println("delta_combined: ", delta_combined)
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

e_msd_1 = percentile(delta_msd, 1)
println("e_msd_1: ", e_msd_1)
e_msd_01 = percentile(delta_msd, 0.1)
println("e_msd_01: ", e_msd_01)

e_combined_1 = percentile(delta_combined, 1)
println("e_combined_1: ", e_combined_1)
e_combined_01 = percentile(delta_combined, 0.1)
println("e_combined_01: ", e_combined_01)

# CALCULATING THE ACCEPTED M' & V' VALUES FOR PLOTTING
accepted_m_si_1 = Float64[]
accepted_m_si_01 = Float64[]
accepted_m_s_1 = Float64[]
accepted_m_s_01 = Float64[]
accepted_m_msd_1 = Float64[]
accepted_m_msd_01 = Float64[]
accepted_m_combined_1 = Float64[]
accepted_m_combined_01 = Float64[]

accepted_v_si_1 = Float64[]
accepted_v_si_01 = Float64[]
accepted_v_s_1 = Float64[]
accepted_v_s_01 = Float64[]
accepted_v_msd_1 = Float64[]
accepted_v_msd_01 = Float64[]
accepted_v_combined_1 = Float64[]
accepted_v_combined_01 = Float64[]

zipped_SI = zip(delta_SI, means, variance)
zipped_S = zip(delta_S, means, variance)
zipped_msd = zip(delta_msd, means, variance)
zipped_combined = zip(delta_combined, means, variance)

# PLOTTING THE POSTERIOR DISTRIBUTION OF THE MEAN STEP LENGTH
# Plot the posterior distribution of the mean step length using S and SI each
# time using 1 and 0.1 percnetiles
# -----------------------------------------------------------------------------
# 1. SI_1
for i in zipped_SI
    if i[1] <= e_SI_1
        push!(accepted_m_si_1, i[2])
        push!(accepted_v_si_1, i[3])
    end
end
x = accepted_m_si_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length", color="salmon")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: PRW: SI_1")

y = accepted_v_si_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: SI_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
sns.jointplot(x="msl", y="variance", data=df, kind="kde")

println("size m': ", size(means))
println("size v': ", size(variance))
println("size accepted_m_si_1: ", size(accepted_m_si_1))
println("size accepted_v_si_1: ", size(accepted_v_si_1))
#------------------------------------------------------------------------------
# # 2. SI_01
# for i in zipped_SI
#     if i[1] <= e_SI_01
#         push!(accepted_m_si_01, i[2])
#         push!(accepted_v_si_01, i[3])
#     end
# end
# x = accepted_m_si_01
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: PRW: SI_0.1")
#
# y = accepted_v_si_01
# fig,ax = PyPlot.subplots()
# sns.distplot(y, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution: PRW: SI_01")
#
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")
#
# println("size m': ", size(means))
# println("size v': ", size(variance))
# println("size accepted_m_si_01: ", size(accepted_m_si_01))
# println("size accepted_v_si_01: ", size(accepted_v_si_01))
# ------------------------------------------------------------------------------
# 3. S_1
for i in zipped_S
    if i[1] <= e_S_1
        push!(accepted_m_s_1, i[2])
        push!(accepted_v_s_1, i[3])
    end
end
x = accepted_m_s_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length", color="salmon")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: PRW: S_1")

y = accepted_v_s_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: S_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
sns.jointplot(x="msl", y="variance", data=df, kind="kde")

println("size m': ", size(means))
println("size v': ", size(variance))
println("size accepted_m_s_1: ", size(accepted_m_s_1))
println("size accepted_v_s_1: ", size(accepted_v_s_1))
# ------------------------------------------------------------------------------
# 4. S_01
# for i in zipped_S
#     if i[1] <= e_S_01
#         push!(accepted_m_s_01, i[2])
#         push!(accepted_v_s_01, i[3])
#     end
# end
# x = accepted_m_s_01
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: PRW: S_01")
#
# y = accepted_v_s_01
# fig,ax = PyPlot.subplots()
# sns.distplot(y, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution: PRW: S_01")
#
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")
#
# println("size m': ", size(means))
# println("size v': ", size(variance))
# println("size accepted_m_s_01: ", size(accepted_m_s_01))
# println("size accepted_v_s_01: ", size(accepted_v_s_01))
# -----------------------------------------------------------------------------
# 5. MSD_1
for i in zipped_msd
    if i[1] <= e_msd_1
        push!(accepted_m_msd_1, i[2])
        push!(accepted_v_msd_1, i[3])
    end
end
x = accepted_m_msd_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length", color="salmon")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: PRW: MSD_1")

y = accepted_v_msd_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: MSD_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
sns.jointplot(x="msl", y="variance", data=df, kind="kde")

println("size m': ", size(means))
println("size v': ", size(variance))
println("size accepted_m_msd_1: ", size(accepted_m_msd_1))
println("size accepted_v_msd_1: ", size(accepted_v_msd_1))
#------------------------------------------------------------------------------
# 6. MSD_01
# for i in zipped_msd
#     if i[1] <= e_msd_01
#         push!(accepted_m_msd_01, i[2])
#         push!(accepted_v_msd_01, i[3])
#     end
# end
# x = accepted_m_msd_01
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: PRW: MSD_0.1_percentile")
#
# y = accepted_v_msd_01
# fig,ax = PyPlot.subplots()
# sns.distplot(y, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution: PRW: MSD_0.1_percentile")
#
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")
#
# println("size m': ", size(means))
# println("size v': ", size(variance))
# println("size accepted_m_msd_01: ", size(accepted_m_msd_01))
# println("size accepted_v_msd_01: ", size(accepted_v_msd_01))
# ------------------------------------------------------------------------------
# 5. COMBINED_1
for i in zipped_combined
    if i[1] <= e_combined_1
        push!(accepted_m_combined_1, i[2])
        push!(accepted_v_combined_1, i[3])
    end
end
x = accepted_m_combined_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length", color="salmon")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: PRW: COMBINED_1")

y = accepted_v_combined_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="Variance")
ax[:set_xlim]([0,2])
ax[:set_title]("Variance Posterior Distribution: PRW: COMBINED_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
sns.jointplot(x="msl", y="variance", data=df, kind="kde")

println("size m': ", size(means))
println("size v': ", size(variance))
println("size accepted_m_combined_1: ", size(accepted_m_combined_1))
println("size accepted_v_combined_1: ", size(accepted_v_combined_1))
# ------------------------------------------------------------------------------
# 6. COMBINED_01
# for i in zipped_combined
#     if i[1] <= e_combined_01
#         push!(accepted_m_combined_01, i[2])
#         push!(accepted_v_combined_01, i[3])
#     end
# end
# x = accepted_m_combined_01
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length", color="salmon")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: PRW: COMBINED_0.1")
#
# y = accepted_v_combined_01
# fig,ax = PyPlot.subplots()
# sns.distplot(y, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution: PRW: COMBINED_0.1")
#
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")
#
# println("size m': ", size(means))
# println("size v': ", size(variance))
# println("size accepted_m_combined_01: ", size(accepted_m_combined_01))
# println("size accepted_v_combined_01: ", size(accepted_v_combined_01))
# ------------------------------------------------------------------------------
# a = sort(accepted_m_si_01)
# fig,ax = PyPlot.subplots()
# accepted_m_100 = a[1:101]
# sns.distplot(a, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: PRW: S_1: Best 100")
#
# b = sort(accepted_v_si_01)
# fig,ax = PyPlot.subplots()
# accepted_v_100 = b[1:101]
# sns.distplot(b, axlabel="Variance")
# ax[:set_xlim]([0,2])
# ax[:set_title]("Variance Posterior Distribution: PRW: S_1: Best 100")
#
# a = x
# b = y
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :variance=>y))
# sns.jointplot(x="msl", y="variance", data=df, kind="kde")




# println("size m': ", size(means))
# println("size accepted_m_si_1: ", size(accepted_m_si_1))
# println("size accepted_m_si_01: ", size(accepted_m_si_01))
# println("size accepted_m_s_1: ", size(accepted_m_s_1))
# println("size accepted_m_s_01: ", size(accepted_m_s_01))
#
# println("size v': ", size(variance))
# println("size accepted_v_si_1: ", size(accepted_v_si_1))
# println("size accepted_v_si_01: ", size(accepted_v_si_01))
# println("size accepted_v_s_1: ", size(accepted_v_s_1))
# println("size accepted_v_s_01: ", size(accepted_v_s_01))
