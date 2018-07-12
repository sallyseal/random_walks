using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

PyPlot.PyObject(PyPlot.axes3D)

# Random walk in 3D
# Using the spherical coordinate system - physics system
# In addition to r and theta we need phi to give us the new point in 3D
# How do we now encorporate bias and persistence into this using the parameters p, b, w?
# Is this even in continuous time or discrete time and how do we change this? Gillespie?

# Initialize vectors with zeros - same number as nsteps
x = zeros(1000)
y = zeros(1000)
z = zeros(1000)

# Initialize a vector the size of the number of steps you want to take
nsteps = zeros(1000)

# Need starting condition for x, y, z as it's 3D - start at the origin
x[1] = 0.0;
y[1] = 0.0;
z[1] = 0.0;

# Iterate through number of steps each time updating x, y and z values and store them in an array to plot
# Start at 2 due to x[i-1]
for i = 2:length(nsteps)

    # need to have a step size so we know how far to travel in the random direction
    # can set a step size or we can randomly sample from some distribution. Sample from normal distribution
    r = rand(Normal())
    # theta is equivalent to latitude and is a value somewhere between -pi and +pi
    # theta = pi*rand() results in clustering of points at poles so we must correct this. It's biased
    # we can't sample uniformally accross theta from 0 - pi we can maybe use inverse cumulitive distribution func.
    # this is the challenge to get unbiased sampling here
    theta = pi*rand()
    # phi is the azimuthal angle in longitude
    phi = 2*pi*rand()

#     println(theta)

    # step size ito x, y, z is
    dx = r*sin(theta)*cos(phi);
    dy = r*sin(theta)*sin(phi);
    dz = r*cos(theta);

    # need to find the new position at the end of this step just taken
    x[i] = x[i-1] + dx
    y[i] = y[i-1] + dy
    z[i] = z[i-1] + dz

end

using PyPlot; const plt = PyPlot

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)

################################################################

# checking the code
# rather than plot random walk just plot set of random points on a sphere - check if we get a sphere
# this code show that there is bias at the poles which is due to the fact that sampling theta is inherantly bias
# generating same number of points at poles as we are at the equator and therefore they cluser more at the poles
# we need to sample theta in another way: inverse cdf method rather than sample uniformally from

# WHY IS THIS NOT WORKING FOR << 1000 POINTS?

using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

y = zeros(1000)
z = zeros(1000)

npoints = zeros(1000)

for i = 2:length(npoints)

    # value of r doesn't matter now
    r = 1
    theta = pi*rand()
    phi = 2*pi*rand()

    x[i] = r*sin(theta)*cos(phi);
    y[i] = r*sin(theta)*sin(phi);
    z[i] = r*cos(theta);

end

# when plotting now we don't want the points connected, but rather just dots to see if we get a sphere
x = x
y = y
z = z

using PyPlot
scatter3D(x, y, z)

######################################################################

using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

# Random walk in 3D - removing the biased sampling of theta at the poles vs. equator

# Initialize vectors with zeros - same number as nsteps
x = zeros(1000)
y = zeros(1000)
z = zeros(1000)

# Initialize a vector the size of the number of steps you want to take
nsteps = zeros(1000)

# Need starting condition for x, y, z as it's 3D - start at the origin
x[1] = 0.0;
y[1] = 0.0;
z[1] = 0.0;

# Iterate through number of steps each time updating x, y and z values and store them in an array to plot
# Start at 2 due to x[i-1]
for i = 2:length(nsteps)

    # need to have a step size so we know how far to travel in the random direction
    # can set a step size or we can randomly sample from some distribution. Sample from normal distribution
    r = rand(Normal())
    # theta is equivalent to latitude and is a value somewhere between -pi and +pi
    # theta = pi*rand() results in clustering of points at poles so we must correct this. It's biased
    # we can't sample uniformally accross theta from 0 - pi we can maybe use inverse cumulitive distribution func.
    # this is the challenge to get unbiased sampling here

    # theta = pi*rand()

    # here is unbiased sampling of theta:
    theta = acos(1-2*rand())

    # phi is the azimuthal angle in longitude
    phi = 2*pi*rand()

#     println(theta)

    # step size ito x, y, z is
    dx = r*sin(theta)*cos(phi);
    dy = r*sin(theta)*sin(phi);
    dz = r*cos(theta);

    # need to find the new position at the end of this step just taken
    x[i] = x[i-1] + dx
    y[i] = y[i-1] + dy
    z[i] = z[i-1] + dz

end

using PyPlot; const plt = PyPlot

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)

#################################################################

# checking the code is no longer biased to theta sampling

# WHY IS THIS NOT WORKING FOR << 1000 POINTS?

using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

y = zeros(1000)
z = zeros(1000)

npoints = zeros(1000)

for i = 2:length(npoints)

    # value of r doesn't matter now
    r = 1
    # arccos = inverse of the cosine function
    theta = acos(1-2*rand())
    phi = 2*pi*rand()

    x[i] = r*sin(theta)*cos(phi);
    y[i] = r*sin(theta)*sin(phi);
    z[i] = r*cos(theta);

end

# when plotting now we don't want the points connected, but rather just dots to see if we get a sphere
x = x
y = y
z = z

using PyPlot
scatter3D(x, y, z)
