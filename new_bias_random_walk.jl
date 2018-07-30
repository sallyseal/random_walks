using Gadfly;
using Distributions;
using PyPlot;
import Plots;

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
  println("source point: ", source)
  println("btheta: ", btheta)
  println("bphi: ", bphi)

  # Perform a RW of nsteps
  for i = 2:length(x)

      # Sample holding time from exponential distribution or another dist?
      t_next_jump = rand(Exponential())
      # Update the time
      t = t+t_next_jump

      # Creating a random point in 3D
      # k = concentration, the higher k, the more biased the random walk
      k = 3
      r = rand(TruncatedNormal(0,1,0,1))
      theta = rand(VonMises(btheta, k),1)      # theta between 0:pi radians
      theta = theta[1]
      phi = rand(VonMises(bphi, k),1)          # phi between 0:2*pi radians
      phi = phi[1]

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
  PyPlot.title("Shape of Random Walk")
