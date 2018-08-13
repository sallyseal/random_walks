using Gadfly;
using Distributions;
using PyPlot;
# using Plotly;
import Plots;

# List of summary statistics to plot
rw_si = Float64[]
rw_si_cart = Float64[]
rw_sinuosity = Float64[]
# rw_msd = Float64[]
# rw_D = Float64[]

# Number of iterations to perform of an nstep random walk
iterations = 1
walkers = zeros(iterations)
for i = 1:length(walkers)

    # Initialize vectors to store the xyz coordinates the size of nsteps
    nsteps = 200
    x = zeros(nsteps)
    y = zeros(nsteps)
    z = zeros(nsteps)

    # Set initial time = 0
    t = 0

    # Create vectors to store variables
    all_r = Float64[]
    all_theta = Float64[]
    all_phi = Float64[]
    time = Float64[]
    holding_time = Float64[]
    all_dtheta = Float64[]
    all_dphi = Float64[]
    all_x = Float64[]
    all_y = Float64[]
    all_z = Float64[]
    turn_angles = Float64[]
    displacements = Float64[]

    # Create starting position of the RW at the origin
    x[1] = 0.0;
    y[1] = 0.0;
    z[1] = 0.0;

    # Sample a random point that will be the source with r = pi
    r = pi
    theta = acos(1-2*rand()) # theta between 0:pi radians
    phi = 2*pi*rand()        # phi between 0:2*pi radians

    source = (r, theta, phi)
    btheta = source[2]
    bphi = source[3]
    # println("source point: ", source)
    # println("btheta: ", btheta)
    # println("bphi: ", bphi)

    # Perform a RW of nsteps
    for i = 2:length(x)

        # Sample holding time from exponential distribution or another dist?
        t_next_jump = rand(Exponential())
        # Update the time
        t = t+t_next_jump

        # Creating a random point in 3D
        # k = concentration, the higher k, the more biased the random walk
        k = 7
        r = rand(TruncatedNormal(0,1,0,1))
        theta = rand(VonMises(btheta, k),1)      # theta between 0:pi radians
        theta = theta[1]
        phi = rand(VonMises(bphi, k),1)          # phi between 0:2*pi radians
        phi = phi[1]
        println(theta)

        # Mapping spherical coordinates onto the cartesian plane
        dx = r*sin(theta)*cos(phi);
        dy = r*sin(theta)*sin(phi);
        dz = r*cos(theta);

        # Updated position
        x[i] = x[i-1] + dx
        y[i] = y[i-1] + dy
        z[i] = z[i-1] + dz

        # Get the current [i] and previous [i-1] coordinates to calculate angle
        # between the 2 vectors = turning angle
        c_1 = x[i], y[i], z[i]
        c_0 = x[i-1], y[i-1], z[i-1]

        # Calculate the turning angle between this vector and previous vector
        turn_angle = acos(vecdot(c_0,c_1)/sqrt(sum(c_1.*c_1)*sum(c_0.*c_0)))

        # Calculate the displacement between i and i-1 and store in vector
        displacement = (x[i-1] - x[i])^2 + (y[i-1] - y[i])^2 + (z[i-1] - z[i])^2
        push!(displacements, displacement)

        # Push to store all values associated with a coordinate
        push!(all_r, r)
        push!(all_theta, theta)
        push!(all_phi, phi)
        push!(time, t)
        push!(holding_time, t_next_jump)
        push!(turn_angles, turn_angle)
        push!(all_x, x[i])
        push!(all_y, y[i])
        push!(all_z, z[i])
    end

    # CALCULATE SUMMARY STATISTICS

    # Straightness Index (D / L):
    # where D = max displacement; L = total path length
    # D = r - r' = sqrt((x-x')^2 + (y-y')^2 + (x-x')^2)
    # SPHERICAL SYSTEM
    theta1 = all_theta[1]
    theta2 = all_theta[end]
    phi1 = all_phi[1]
    phi2 = all_phi[end]
    r1 = all_r[1]
    r2 = all_r[end]

    D = r1^2 + r2^2 -
        2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1 - phi2) + cos(theta1)*cos(theta2))
    d = sqrt(D)
    L = sum(all_r)
    si = d/L

    # CARTESIAN SYSTEM
    x1 = all_x[1]
    x2 = all_x[end]
    y1 = all_y[1]
    y2 = all_y[end]
    z1 = all_z[1]
    z2 = all_z[end]

    disp = (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
    disp_sqrt = sqrt(disp)
    si_cart = disp_sqrt/L

    # Sinuosity Index: measures path deviation locally s prop sd/mur
    # where sd = standard dev of turn angle distribution
    # mur = mean step length
    mur = mean(all_r)
    sd = std(turn_angles[2:end])
    sinuosity = sd/mur

    # Push values to a list to store them for later statistics
    push!(rw_si, si)
    push!(rw_si_cart, si_cart)
    push!(rw_sinuosity, sinuosity)

    # PLOTTING

    # Plotting RW for each iteration
    # Uncomment the below if you want to visualise each walk
    using PyPlot; const plt = PyPlot
    PyPlot.PyObject(PyPlot.axes3D)

    x = x
    y = y
    z = z

    fig = plt.figure()
    ax = fig[:add_subplot](111, projection="3d")
    ax[:plot](x, y, z)
    # PyPlot.title("Shape of Random Walk")
    PyPlot.xlabel("x")
    PyPlot.ylabel("y")
    PyPlot.zlabel("z")
end

# Calculate the mean of summary statistics
rw_si_mu = mean(rw_si)
rw_si_cart_mu = mean(rw_si_cart)
rw_sinuosity_mu = mean(rw_sinuosity)
println("rw straightness index average: ", rw_si_mu)
println("rw_cartesian straightness index average: ", rw_si_cart_mu)
println("rw sinuosity average: ", rw_sinuosity_mu)
# println("rw sinuosity: ", rw_sinuosity)

# Plotting distributions of straightness index
# a = rw_si
# plot1 = PyPlot.plt[:hist](a)
# PyPlot.xlabel("Straightness Index")
# PyPlot.title("Random Walk Straightness Index Histogram")

# a_cart = rw_si_cart
# plot1 = PyPlot.plt[:hist](a_cart)
# PyPlot.xlabel("Straightness Index")
# PyPlot.title("Randon Walk Straightness Index Cartesian")

# Plotting distributions of the sinuosity
# b = rw_sinuosity
# plot2 = PyPlot.plt[:hist](b)
# PyPlot.xlabel("Sinuosity")
# PyPlot.title("Randon Walk Sinuosity Histogram")
