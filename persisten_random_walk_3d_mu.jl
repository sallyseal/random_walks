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

# Initialize vectors
x = zeros(1000)
y = zeros(1000)
z = zeros(1000)

# Set initial time = 0 and have a total
t = 0
total_time = 1000

# Create vectors to store r, theta, phi, time, holding time for each xyz coordinate
all_r = Float64[]
all_theta = Float64[]
all_phi = Float64[]
time = Float64[]
holding_time = Float64[]

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
r = rand(TruncatedNormal(0,1,0,1))
theta = acos(1-2*rand()) # theta between 0:pi radians
phi = 2*pi*rand()        # phi between 0:2*pi radians

# FOR THE PERSISTENCE: variance
sigma_t = 0.4 # Can control the tightness/spread of the distribution by altering
sigma_p = 0.4 # Can control the tightness/spread of the distribution by altering

# Perform simulation while t is <= total time of the reaction
while t <= total_time

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
        dist_theta = TruncatedNormal(theta, sigma_t, lower_t, upper_t)
        dist_phi = TruncatedNormal(phi, sigma_p, lower_p, upper_p)

        # Randomly sample from the distributions to get updated theta and phi to
        # create next point in 3D space
        theta = rand(dist_theta)
        phi = rand(dist_phi)
        r = rand(TruncatedNormal(0,1,0,1))

        # Calculate dtheta and dphi: angle between new and old theta and phi
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

        # Push to store all values associated with a coordinate
        push!(all_r, r)
        push!(all_theta, theta)
        push!(all_phi, phi)
        push!(time, t)
        push!(holding_time, t_next_jump)
        push!(all_dtheta, dtheta)
        push!(all_dphi, dphi)
    end
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

norm = r1^2 + r2^2 -
    2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1 - phi2) + cos(theta1)*cos(theta2))
l = sum(all_r)
s_index = norm/l

println(s_index)

# Sinuosity Index: measures path deviation locally s prop sigma theta/ mux
# where sigma theta = standard dev of turn angle distribution
# mux is is mean step length
mux = mean(all_r)
sigma_theta = std(all_dtheta)
sigma_phi = std(all_dphi)
println(sigma_theta)

# Plotting
using PyPlot; const plt = PyPlot
PyPlot.PyObject(PyPlot.axes3D)

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)
