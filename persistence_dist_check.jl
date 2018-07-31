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

# Create starting position at the origin
x[1] = 0.0;
y[1] = 0.0;
z[1] = 0.0;

# Sample first random point in 3D
r = 1
theta = acos(1-2*rand()) # theta between 0:pi radians
phi = 2*pi*rand()        # phi between 0:2*pi radians

# FOR THE PERSISTENCE
mu_t = 0
sigma_t = 0.4 # Can control the tightness/spread of the distribution by altering
mu_p = 0
sigma_p = 0.4 # Can control the tightness/spread of the distribution by altering

# Create the distributions for theta and phi - change with which to update theta and phi
dist_theta = TruncatedNormal(mu_t, sigma_t, -0.5*pi, 0.5*pi)
dist_phi = TruncatedNormal(mu_p, sigma_p, -pi, pi)

# Perform simulation while t is <= total time of the reaction
while t <= total_time

    for i = 2:length(x)

        # Sample holding time from exponential distribution or another dist?
        t_next_jump = rand(Exponential())
        # Update the time
        t = t+t_next_jump

        # Randomly sample from the distributions to get the variance of the next angle
        dtheta = rand(dist_theta)
        dphi = rand(dist_phi)

        # Update theta and phi with variance sampled from distributions above
        theta += dtheta
        phi += dphi
        r = 1

        x[i] = r*sin(theta)*cos(phi);
        y[i] = r*sin(theta)*sin(phi);
        z[i] = r*cos(theta);

        # Push to store all values associated with a coordinate
        push!(all_r, r)
        push!(all_theta, theta)
        push!(all_phi, phi)
        push!(time, t)
        push!(holding_time, t_next_jump)

    end
end

# when plotting now we don't want the points connected, but rather just dots to see if we get a sphere
x = x
y = y
z = z

using PyPlot
scatter3D(x, y, z)
