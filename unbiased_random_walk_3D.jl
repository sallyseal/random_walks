using Gadfly;
using Distributions;
using PyPlot;
using Plotly;

# List of summary statistics to plot
rw_si = Float64[]
rw_sinuosity = Float64[]
# rw_msd = Float64[]
# rw_D = Float64[]

# Number of iterations to perform of an nstep random walk
iterations = 5
walkers = zeros(iterations)
for i = 1:length(walkers)

    # Initialize vectors to store the xyz coordinates the size of nsteps
    nsteps = 100
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
    # all_msd = Float64[]

    # Create starting position of the RW at the origin
    x[1] = 0.0;
    y[1] = 0.0;
    z[1] = 0.0;

    # Perform a RW of nsteps
    for i = 2:length(x)

        # Sample holding time from exponential distribution or another dist?
        t_next_jump = rand(Exponential())
        # Update the time
        t = t+t_next_jump

        # Creating a random point in 3D
        # Same distribution used in Cell Press Paper
        r = rand(TruncatedNormal(0,1,0,1))
        theta = acos(1-2*rand()) # theta between 0:pi radians
        phi = 2*pi*rand()        # phi between 0:2*pi radians

        # msd = r[i]^2 + r[i-1]^2 - 2*r[i]*r[i-1]*(sin(theta[i-1])*sin(theta[i])*cos(phi[i-1] - phi[i]) + cos(theta[i-1])*cos(theta[i]))
        # msd = msd^2

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

        # Push to store all values associated with a coordinate
        push!(all_r, r)
        push!(all_theta, theta)
        push!(all_phi, phi)
        push!(time, t)
        push!(holding_time, t_next_jump)
        push!(turn_angles, turn_angle)
        # push!(all_msd, msd)
    end

    # CALCULATE SUMMARY STATISTICS

    # Straightness Index (D / L):
    # where D = max displacement; L = total path length
    # D = r - r' = sqrt((x-x')^2 + (y-y')^2 + (x-x')^2)
    theta1 = all_theta[1]
    theta2 = all_theta[end]
    phi1 = all_phi[1]
    phi2 = all_phi[end]
    r1 = all_r[1]
    r2 = all_r[end]

    D = r1^2 + r2^2 -
        2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1 - phi2) + cos(theta1)*cos(theta2))
    L = sum(all_r)
    si = D/L
    # println("straightness index:", si)

    # Sinuosity Index: measures path deviation locally s prop sd/mur
    # where sd = standard dev of turn angle distribution
    # mur = mean step length
    mur = mean(all_r)
    sd = std(turn_angles[2:end])
    sinuosity = sd/mur
    # println("Sinuosity: ", sinuosity)

    push!(rw_si, si)
    push!(rw_sinuosity, sinuosity)
    println("rw straightness indexes: ", rw_si)
    println("rw sinuosity: ", rw_sinuosity)


    # Plotting RW uncomment the below if you want to visualise each walk
    # using PyPlot; const plt = PyPlot
    # PyPlot.PyObject(PyPlot.axes3D)
    #
    # x = x
    # y = y
    # z = z
    #
    # fig = plt.figure()
    # ax = fig[:add_subplot](111, projection="3d")
    # ax[:plot](x, y, z)

    # Plotting distributions of straightness index & sinuosity
    a = rw_si
    plot1 = plt[:hist](a)

    b = rw_sinuosity
    plot2 = plt[:hist](b)


    # Plotting msd vs time
    # Maybe do log of this/ work out d and use this as summary stat
end
