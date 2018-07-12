using Gadfly;
using Distributions;
using PyPlot;
using Plotly;
PyPlot.PyObject(PyPlot.axes3D)

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
