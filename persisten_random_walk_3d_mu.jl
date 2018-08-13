# This script creates a persistent random walk where the new theta and phi are sampled
# by creating a truncated normal distribution and updating mu with the old theta
# and phi values, respectively. Each time the bounds are also updated accordingly
# so that theta is between 0:pi and phi 0:2pi

# We create 2 random points (with the second points theta and phi relying on the
# first point, but dtheta and dphi are calculated retrospectively) i.e. we don't
# calculate a dtheta and dphi and then apply it to create a new point
using Gadfly;
using Distributions;
using PyPlot;
using Plotly;
using StatPlots;

# List of summary statistics to plot
prw_si = Float64[]
prw_si_cart = Float64[]
prw_sinuosity = Float64[]

iterations = 1
walkers = zeros(iterations)
for i = 1:length(walkers)

    # Initialize vectors the size of nsteps
    nsteps = 150
    x = zeros(nsteps)
    y = zeros(nsteps)
    z = zeros(nsteps)

    # Set initial time = 0
    t = 0

    # Create vectors to store r, theta, phi, time, holding time for each xyz coordinate
    all_r = Float64[]
    all_theta = Float64[]
    all_phi = Float64[]
    all_x = Float64[]
    all_y = Float64[]
    all_z = Float64[]
    time = Float64[]
    holding_time = Float64[]
    turn_angles = Float64[]

    # Bounds for distributions
    lower_t = 0
    upper_t = pi
    lower_p = 0
    upper_p = 2*pi

    # Persistence Angle
    all_dtheta = Float64[]
    all_dphi = Float64[]

    # Create starting position at the origin
    x[1] = 0.0;
    y[1] = 0.0;
    z[1] = 0.0;

    # Sample first random point in 3D
    r = rand(TruncatedNormal(0,1,0,1)) # Adam uses log normal?
    theta = acos(1-2*rand()) # theta between 0:pi radians
    phi = 2*pi*rand()        # phi between 0:2*pi radians

    # FOR THE PERSISTENCE: variance
    sigma_t = 0.2 # Can control the tightness/spread of the distribution by altering
    sigma_p = 0.2 # Can control the tightness/spread of the distribution by altering

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
        dist_theta = TruncatedNormal(theta, sigma_t, lower_t, upper_t)
        dist_phi = TruncatedNormal(phi, sigma_p, lower_p, upper_p)

        # Randomly sample from the distributions to get updated theta and phi to
        # create next point in 3D space
        theta = rand(dist_theta)
        phi = rand(dist_phi)
        r = rand(TruncatedNormal(0,1,0,1))

        # Calculate dtheta and dphi: angle between new and old theta and phi
        # Or should this rahter be dot product
        dtheta = mu_t - theta
        dphi = mu_p - phi

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
        push!(all_r, r)
        push!(all_theta, theta)
        push!(all_phi, phi)
        push!(time, t)
        push!(holding_time, t_next_jump)
        push!(all_dtheta, dtheta)
        push!(all_dphi, dphi)
        push!(all_x, x[i])
        push!(all_y, y[i])
        push!(all_z, z[i])
        push!(turn_angles, turn_angle)
    end

    # CALCULATE SUMMARY STATISTICS

    # Straightness Index: D/L where D = max displacement; L = total path length
    # D: r - r' = sqrt((x-x')^2 + (y-y')^2 + (x-x')^2)
    theta1 = all_theta[1]
    theta2 = all_theta[end]
    phi1 = all_phi[1]
    phi2 = all_phi[end]
    r1 = all_r[1]
    r2 = all_r[end]

    D = r1^2 + r2^2 -
        2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1 - phi2) + cos(theta1)*cos(theta2))
    d = sqrt(D)
    L = sum(all_r)
    si = d/L

    # CARTESIAN SYSTEM
    x1 = all_x[1]
    x2 = all_x[end]
    y1 = all_y[1]
    y2 = all_y[end]
    z1 = all_z[1]
    z2 = all_z[end]

    disp = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
    disp_sqrt = sqrt(disp)
    si_cart = disp_sqrt/L

    # Sinuosity Index: measures path deviation locally s prop sd/ mur
    # where sd = standard dev of turn angle distribution
    # mur is is mean step length
    mur = mean(all_r)
    sd = std(turn_angles[2:end])
    sinuosity = sd/mur

    # Push values to a list to store them for later statistics
    push!(prw_si, si)
    push!(prw_si_cart, si_cart)
    push!(prw_sinuosity, sinuosity)

    # Plotting each RW - uncomment if want to see this
    using PyPlot; const plt = PyPlot
    PyPlot.PyObject(PyPlot.axes3D)

    x = x
    y = y
    z = z

    fig = plt.figure()
    ax = fig[:add_subplot](111, projection="3d")
    ax[:plot](x, y, z)
    # PyPlot.title("Shape of Persistent Random Walk")
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    PyPlot.zlabel("z")
end

# Calculate the mean of summary statistics
prw_si_mu = mean(prw_si)
prw_si_cart_mu = mean(prw_si_cart)
prw_sinuosity_mu = mean(prw_sinuosity)
println("prw straightness index average: ", prw_si_mu)
println("prw_cartesian straightness index average: ", prw_si_cart_mu)
println("prw sinuosity average: ", prw_sinuosity_mu)

# Plotting distributions of straightness index
# a = prw_si
# plot1 = PyPlot.plt[:hist](a)
# PyPlot.xlabel("Straightness Index")
# PyPlot.title("Persistent Randon Walk Straightness Index Histogram")

# a_cart = prw_si_cart
# plot1 = PyPlot.plt[:hist](a_cart)
# PyPlot.xlabel("Straightness Index")
# PyPlot.title("Randon Walk Straightness Index Cartesian")

# Plotting distributions of the sinuosity
# b = prw_sinuosity
# plot2 = PyPlot.plt[:hist](b)
# PyPlot.xlabel("Sinuosity")
# PyPlot.title("Persistent Randon Walk Sinuosity Histogram")
