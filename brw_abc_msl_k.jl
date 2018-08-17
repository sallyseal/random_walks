# This script does the same as script rw_simulations.jl
# Here we attempt to infer k of the theta an phi von Mises distributions
# This is a way to infer how bias the random walk is
# Sample k from a uniform distribution between 0 and 20?
# Check when the rw looks ver biased and then use that k value as the max (20)
# Keep all other parameters constant between mock and simulated data

using Distributions;
using PyPlot;
using StatsBase;
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport seaborn as sns
using PyCall, PyPlot; @pyimport pandas as pd

# Generate the mock data (10x RWs of 100 steps each) and get summary statistics
########## MOCK DATA ##########

# Create vectors to store the average SI and S for 10x RWs
SI_av = Float64[]
S_av = Float64[]

# Set the level or bias higher by increasing k
k = 10
msl = 0.5

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

    # Sample a random point that will be the source with r = pi
    r = pi
    theta = acos(1-2*rand()) # theta between 0:pi radians
    phi = 2*pi*rand()        # phi between 0:2*pi radians

    source = (r, theta, phi)
    btheta = source[2]
    bphi = source[3]

    # Perform a RW of nsteps
    for i = 2:length(x)

        # Sample holding time from exponential distribution or another dist?
        t_next_jump = rand(Exponential())
        # Update the time
        t = t+t_next_jump

        # Creating a random point in 3D
        # k = concentration, the higher k, the more biased the random walk
        # k = 1
        # Here we want to try and infer the mean step length
        r = rand(TruncatedNormal(msl, 0.1, 0, 1))
        theta = rand(VonMises(btheta, k),1)      # theta between 0:pi radians
        theta = theta[1]
        phi = rand(VonMises(bphi, k),1)          # phi between 0:2*pi radians
        phi = phi[1]

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
ks = Float64[]

# Repeat simulation 10 000x
for i in 1:100000

    # Generate the simulated data (10x RWs of 100 steps each) and get summary stats
    ########## SIMULATED DATA ##########

    # Create vectors to store the average SI and S for 10x RWs
    SI_prime_av = Float64[]
    S_prime_av = Float64[]

    # Sample k from uniform dist between 0 & 20? save value to means
    # Check how to use Truncated uniform distribution
    m = rand()
    push!(means, m)
    k_prime = rand(Uniform(0,10))
    push!(ks, k_prime)

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

        # Sample a random point that will be the source with r = pi
        r = pi
        theta = acos(1-2*rand()) # theta between 0:pi radians
        phi = 2*pi*rand()        # phi between 0:2*pi radians

        source = (r, theta, phi)
        btheta = source[2]
        bphi = source[3]

        # Perform a RW of nsteps
        for i = 2:length(x)

            # Sample holding time from exponential distribution or another dist?
            t_next_jump = rand(Exponential())
            # Update the time
            t = t+t_next_jump

            # Creating a random point in 3D
            # k = concentration, the higher k, the more biased the random walk
            # k = 1
            # Here we want to try and infer the mean step length
            r = rand(TruncatedNormal(m, 0.1, 0, 1))
            theta = rand(VonMises(btheta, k_prime),1)      # theta between 0:pi radians
            theta = theta[1]
            phi = rand(VonMises(bphi, k_prime),1)          # phi between 0:2*pi radians
            phi = phi[1]

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
    SI_prime_av = mean(SI_prime_av)
    S_prime_av = mean(S_prime_av)

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

# CALCULATING THE ACCEPTED M'& K' VALUES FOR PLOTTING
accepted_m_si_1 = Float64[]
accepted_m_si_01 = Float64[]
accepted_m_s_1 = Float64[]
accepted_m_s_01 = Float64[]

accepted_k_si_1 = Float64[]
accepted_k_si_01 = Float64[]
accepted_k_s_1 = Float64[]
accepted_k_s_01 = Float64[]

zipped_SI = zip(delta_SI, means, ks)
zipped_S = zip(delta_S, means, ks)

# PLOTTING THE POSTERIOR DISTRIBUTION OF THE MEAN STEP LENGTH
# Plot the posterior distribution of the mean step length using S and SI each
# time using 1 and 0.1 percnetiles
# -----------------------------------------------------------------------------
# 1. SI_1
for i in zipped_SI
    if i[1] <= e_SI_1
        push!(accepted_m_si_1, i[2])
        push!(accepted_k_si_1, i[3])
    end
end
x = accepted_m_si_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: BRW: SI_1")

y = accepted_k_si_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="K")
ax[:set_xlim]([0,10])
ax[:set_title]("K Posterior Distribution: BRW: SI_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :k=>y))
sns.jointplot(x="msl", y="k", data=df, kind="kde")

println("size m': ", size(means))
println("size k': ", size(ks))
println("size accepted_m_si_1: ", size(accepted_m_si_1))
println("size accepted_k_si_1: ", size(accepted_k_si_1))
#------------------------------------------------------------------------------
#
# # 2. SI_01
for i in zipped_SI
    if i[1] <= e_SI_01
        push!(accepted_m_si_01, i[2])
        push!(accepted_k_si_01, i[3])
    end
end
x = accepted_m_si_01
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: BRW: SI_0.1")

y = accepted_k_si_01
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="K")
ax[:set_xlim]([0,10])
ax[:set_title]("K Posterior Distribution: BRW: SI_01")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :k=>y))
sns.jointplot(x="msl", y="k", data=df, kind="kde")

println("size m': ", size(means))
println("size k': ", size(ks))
println("size accepted_m_si_01: ", size(accepted_m_si_01))
println("size accepted_k_si_01: ", size(accepted_k_si_01))
# ------------------------------------------------------------------------------
#
# 3. S_1
for i in zipped_S
    if i[1] <= e_S_1
        push!(accepted_m_s_1, i[2])
        push!(accepted_k_s_1, i[3])
    end
end
x = accepted_m_s_1
fig,ax = PyPlot.subplots()
sns.distplot(x, axlabel="Mean Step Length")
ax[:set_xlim]([0,1])
ax[:set_title]("Mean Step Length Posterior Distribution: BRW: S_1")

y = accepted_k_s_1
fig,ax = PyPlot.subplots()
sns.distplot(y, axlabel="K")
ax[:set_xlim]([0,10])
ax[:set_title]("K Posterior Distribution: BRW: S_1")

fig,ax = PyPlot.subplots()
df = pd.DataFrame(data=Dict(:msl=>x, :k=>y))
sns.jointplot(x="msl", y="k", data=df, kind="kde")

println("size m': ", size(means))
println("size k': ", size(ks))
println("size accepted_m_s_1: ", size(accepted_m_s_1))
println("size accepted_k_s_1: ", size(accepted_k_s_1))
# ------------------------------------------------------------------------------

# 4. S_01
# for i in zipped_S
#     if i[1] <= e_S_01
#         push!(accepted_m_s_01, i[2])
#         push!(accepted_k_s_01, i[3])
#     end
# end
# x = accepted_m_s_01
# fig,ax = PyPlot.subplots()
# sns.distplot(x, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: BRW: S_01")
#
# y = accepted_k_s_01
# fig,ax = PyPlot.subplots()
# sns.distplot(y, axlabel="K")
# ax[:set_xlim]([0,10])
# ax[:set_title]("K Posterior Distribution: BRW: S_01")
#
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :k=>y))
# sns.jointplot(x="msl", y="k", data=df, kind="kde")
#
# println("size m': ", size(means))
# println("size k': ", size(ks))
# println("size accepted_m_s_01: ", size(accepted_m_s_01))
# println("size accepted_k_s_01: ", size(accepted_k_s_01))
# -----------------------------------------------------------------------------

# Plot for plotting the top 100 or top 1000 particles in a particular run
# a = sort(accepted_m_si_1)
# fig,ax = PyPlot.subplots()
# accepted_m_100 = a[1:101]
# sns.distplot(a, axlabel="Mean Step Length")
# ax[:set_xlim]([0,1])
# ax[:set_title]("Mean Step Length Posterior Distribution: BRW: Si_1: Best 100")
#
# b = sort(accepted_k_si_1)
# fig,ax = PyPlot.subplots()
# accepted_k_100 = b[1:101]
# sns.distplot(b, axlabel="K")
# ax[:set_xlim]([0,2])
# ax[:set_title]("K Posterior Distribution: BRW: Si_1: Best 100")
#
# a = x
# b = y
# fig,ax = PyPlot.subplots()
# df = pd.DataFrame(data=Dict(:msl=>x, :k=>y))
# sns.jointplot(x="msl", y="k", data=df, kind="kde")
# ------------------------------------------------------------------------------

# println("size m': ", size(means))
# println("size accepted_m_si_1: ", size(accepted_m_si_1))
# println("size accepted_m_si_01: ", size(accepted_m_si_01))
# println("size accepted_m_s_1: ", size(accepted_m_s_1))
# println("size accepted_m_s_01: ", size(accepted_m_s_01))
#
# println("size k': ", size(ks))
# println("size accepted_k_si_1: ", size(accepted_k_si_1))
# println("size accepted_k_si_01: ", size(accepted_k_si_01))
# println("size accepted_k_s_1: ", size(accepted_k_s_1))
# println("size accepted_k_s_01: ", size(accepted_k_s_01))
