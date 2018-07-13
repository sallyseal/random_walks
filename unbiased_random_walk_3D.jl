using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

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

# Perform simulation while t is <= total time of the reaction
while t <= total_time

    for i = 2:length(x)

        # Sample holding time from exponential distribution or another dist?
        t_next_jump = rand(Exponential())
        # Update the time
        t = t+t_next_jump

        # Creating a random point in 3D
        # Step size = r
        # Same distribution used in Cell Press Paper
        r = rand(TruncatedNormal(0,1,0,1))
        theta = acos(1-2*rand()) # theta between 0:pi radians # confirm
        phi = 2*pi*rand()        # phi between 0:2*pi radians

        # mapping spherical coordinates onto the cartesian plane
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

# Plotting
using PyPlot; const plt = PyPlot
PyPlot.PyObject(PyPlot.axes3D)

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)
