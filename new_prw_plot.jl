using Gadfly;
using Distributions;
using PyPlot;
# using Plotly;
using StatPlots;

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

# Create starting position at the origin
x[1] = 0.0;
y[1] = 0.0;
z[1] = 0.0;

# Sample first random point in 3D
r = rand(TruncatedNormal(0,1,0,1))
theta = acos(1-2*rand()) # theta between 0:pi radians
phi = 2*pi*rand()        # phi between 0:2*pi radians

# This is what we want to try to INFER
# FOR THE PERSISTENCE: variance of theta and phi distributions
k = 20 # Can control the tightness/spread of the distribution by altering
msl = 0.5

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
    dist_theta = VonMises(theta, k)
    dist_phi = VonMises(phi, k)

    # Randomly sample from the distributions to get updated theta and phi to
    # create next point in 3D space
    theta = rand(dist_theta,1)
    theta = theta[1]
    phi = rand(dist_phi,1)
    phi = phi[1]
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

# Plotting
using PyPlot; const plt = PyPlot
PyPlot.PyObject(PyPlot.axes3D)

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)
