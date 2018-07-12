# Now try sampling from exponential distribution for step size/ cell displacement velocity

using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

PyPlot.PyObject(PyPlot.axes3D)

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
    r = rand(Exponential())
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
