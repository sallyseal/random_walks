# Cleaned code:

using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

PyPlot.PyObject(PyPlot.axes3D)

x = zeros(1000)
y = zeros(1000)
z = zeros(1000)

# Number of steps in the RW
nsteps = zeros(1000)

# Need starting condition - start at the origin
x[1] = 0.0;
y[1] = 0.0;
z[1] = 0.0;

# Iterate through number of steps each time updating x, y and z values
for i = 2:length(nsteps)

    # Should this be truncated at 0.1?
    r = rand(Exponential())
    theta = acos(1-2*rand())
    phi = 2*pi*rand()

    # step size ito x, y, z is
    dx = r*sin(theta)*cos(phi);
    dy = r*sin(theta)*sin(phi);
    dz = r*cos(theta);

    # Updated position at the end of the step just taken
    x[i] = x[i-1] + dx
    y[i] = y[i-1] + dy
    z[i] = z[i-1] + dzfgvc

end

# Plotting
using PyPlot; const plt = PyPlot

x = x
y = y
z = z

fig = plt.figure()
ax = fig[:add_subplot](111, projection="3d")
ax[:plot](x, y, z)

# to impliment bias and persistence: function is currently sample random number. Replace this with sample from
# wrapped normal distribution (Julia doesn't have this distribution built in - von misses?)

# for bias:
# mean = beta = angle of bias = vector pointing from cell to attractant
# variance = sigma = -2log(b) where b is somewhere between 0 and 1 (can we simply set b to something between 0-1?)

# for persistence:
# mean = angle direction of previous step
# variance = sigma = -2log(p)

# probability of bias = w
# probability of persistence = 1-w
