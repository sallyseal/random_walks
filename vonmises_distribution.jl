# checking the code is no longer biased to theta sampling

# WHY IS THIS NOT WORKING FOR << 1000 POINTS?

using Gadfly;
using Distributions;
using PyPlot;
using Plotly;
PyPlot.PyObject(PyPlot.axes3D)

nsteps = 1000

x = zeros(nsteps)
y = zeros(nsteps)
z = zeros(nsteps)

npoints = zeros(1000)

for i = 2:length(x)

    k = 5
    # value of r doesn't matter now
    r = 1
    # arccos = inverse of the cosine function
    theta = rand(VonMises(1, k),1)
    phi = rand(VonMises(1, k),1)
    theta = theta[1]
    phi = phi[1]

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
