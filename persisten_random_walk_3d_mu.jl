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
sigma_t = 0.1 # Can control the tightness/spread of the distribution by altering
sigma_p = 0.1 # Can control the tightness/spread of the distribution by altering

# Perform simulation while t is <= total time of the reaction


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
# Only theta from wrapped normal distribution right?
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
# println("turn_angle: ", turn_angle)

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

# println(size(all_dphi))
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

println("Straightness Index PRW: ", s_index)

# Sinuosity Index: measures path deviation locally s prop sd/ mux
# where sd = standard dev of turn angle distribution
# mux is is mean step length
mux = mean(all_r)
sd = std(turn_angles[2:end])
println("sd_turn_angle: ", sd)
sinuosity = sd/mux
println("Sinuosity PRW: ", sinuosity)

# Plotting
using PyPlot; const plt = PyPlot
PyPlot.PyObject(PyPlot.axes3D)

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)
