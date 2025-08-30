using LinearAlgebra
using DelimitedFiles
using BenchmarkTools
using DifferentialEquations


function read_txt_file(file_path::AbstractString)
    # Initialize an empty array to store the data
    data_array = []

    # Open the file and read its contents
    open(file_path, "r") do file
        for line in eachline(file)
            # Split each line into individual values
            values = split(line)

            # Convert the values to Float64 and append to the array
            push!(data_array, [parse(Float64, values[1]), parse(Float64, values[2])])
        end
    end

    # Convert the array to a Julia array
    data_array = Array(data_array)

    return data_array
end


## --------------------------------------------------------------- FUNCTIONS FOR RK IMPLEMENTATION --------------------------------------------------------------- ##

## this function calculates r_{ij}^2 ##
function r_squared(pos_Real, i, j)

    return (pos_Real[i,1] - pos_Real[j,1])^2 + (pos_Real[i,2] - pos_Real[j,2])^2

end


## this function calculates and returns the vector that stores the x- and y-coordinates of the ith image vortex - indexing this output by 1 gives teh x-coordinate while index 2 gives the y-coordinate ##
function posImage(pos_Real, i, R)

    return (R^2)*pos_Real[i,:]/(pos_Real[i,1]^2 + pos_Real[i,2]^2)

end


## this function returns and calculates r_bar_squared for the ith vortex ##
function r_bar_squared(pos_Real, i, j)

    return (pos_Real[i,1] - posImage(pos_Real,j,R)[1])^2 + (pos_Real[i,2] - posImage(pos_Real,j,R)[2])^2

end


## this function returns and calculates the right-hand side of the ODE for x and y (with damping) in an array with each row representing a different point vortex and the first and second columns representing the x- and y-coordinates, respectively ##
function derivative(pos_Real, gamma, damping_Coeff, Omega)

    (n, m) = size(pos_Real);
    superfluidVelocity = zeros(Float64, n, m);
    vortexVelocity = zeros(Float64, n, m);

    # this part adds the term that accounts for the rotating frame to the local superfluid velocity
    superfluidVelocity[:,1] = superfluidVelocity[:,1] - Omega*pos_Real[:,2];
    superfluidVelocity[:,2] = superfluidVelocity[:,2] + Omega*pos_Real[:,1];

    # this part will first calculate the local superfluid velocity, v_i, for each point vortex location, and then will add the correction for damping
    for i in 1:n
        ## THE COMMENTED OUT STUFF IN THIS LOOP IS FOR THE BOUNDARY CONDITIONS, WHICH WE WILL IGNORE FOR THE SIMULATION
        for j in 1:n
            if j==i
                continue
               # vec[i,1] = vec[i,1] + (1/2*pi)*gamma[j]*(pos_Real[i,2] - posImage(pos_Real,j,R)[2])/r_bar_squared(pos_Real,i,j) # calculates the x-component
               # vec[i,2] = vec[i,2] - (1/2*pi)*gamma[j]*(pos_Real[i,1] - posImage(pos_Real,j,R)[1])/r_bar_squared(pos_Real,i,j) # calculates the y-component
            else
                superfluidVelocity[i,1] = superfluidVelocity[i,1] + (1/2*pi)*gamma*(pos_Real[i,2] - pos_Real[j,2])/r_squared(pos_Real,i,j); # + (1/2*pi)*gamma[j]*(pos_Real[i,2] - posImage(pos_Real,j,R)[2])/r_bar_squared(pos_Real,i,j) # calculates the x-component
                superfluidVelocity[i,2] = superfluidVelocity[i,2] - (1/2*pi)*gamma*(pos_Real[i,1] - pos_Real[j,1])/r_squared(pos_Real,i,j); # - (1/2*pi)*gamma[j]*(pos_Real[i,1] - posImage(pos_Real,j,R)[1])/r_bar_squared(pos_Real,i,j) # calculates the x-component
            end
        end
        
    end

    # this part calculates the equations of motion for the actual vortices from the local superfluid velocity calculated above
    vortexVelocity[:,1] = superfluidVelocity[:,1] - damping_Coeff*superfluidVelocity[:,2];
    vortexVelocity[:,2] = superfluidVelocity[:,2] + damping_Coeff*superfluidVelocity[:,1];

    return vortexVelocity

end


## this function does the integration step for Runge-Kutta - calculates the part we need to add to each poisition to calculate the positions at the next time step (this will return an array, with each row for a different point vortex, and then the two columns represent x and y) ##
function RK(pos_Real, gamma, damping_Coeff, Omega, dt)
    
    # calculates f0 as just the right-hand side of the ODE
    f0 = derivative(pos_Real, gamma, damping_Coeff, Omega);

    # calculates f1, which is the right-hand side of the ODE at the current position plus 0.5*dt*f0
    f1 = derivative(pos_Real + 0.5*dt*f0, gamma, damping_Coeff, Omega);

    # calculates f2, which is the right-hand side of the ODE at the current position plus 0.5*dt*f1
    f2 = derivative(pos_Real + 0.5*dt*f1, gamma, damping_Coeff, Omega);

    # calculates f3, which is the right-hand side of the ODE at the current position plus 0.5*dt*f1
    f3 = derivative(pos_Real + dt*f2, gamma, damping_Coeff, Omega);

    return (dt/6)*(f0 + 2*f1 + 2*f2 + f3)

end


## --------------------------------------------------------------- DEFINING CONSTANTS AND INITIALISING VARIABLES --------------------------------------------------------------- ##


#= global positions = read_txt_file(raw"C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Code\triangular_lattice.txt"); # triangular lattice initial configuration
NumVortices = size(positions,1);
global vortexPosReal = zeros(Float64, NumVortices, 2);

for j=1:NumVortices
    vortexPosReal[j,1] = positions[j][1];
    vortexPosReal[j,2] = positions[j][2];
end

println(NumVortices) =#



NumVortices = 1000; # this is the number of vortices in the simulation
dt = 0.001; # time step
tmax = 10; # maximum integration time
NumTime = floor(Int,tmax/dt); # this is the number of time steps

R = 1.0; # radius of superfluid

# gamma = ones(Float64, (NumVortices, 1)); # circulation of each point vortex is generated as being equal to 1 - this is from when I was being general and allowing the circulation of each vortex to be different
gamma = 1; # sets the circulation of each vortex to be equal to 1 - just a scalar, rather than what I 

DampingCoeff = 0.3; # damping coefficient
Omega = 5; # angular velocity of the rotating frame

global vortexPosReal = zeros(Float64, (NumVortices, 2)) # initialises array to store the initial point vortex positions

## ----------------------------------------------------------- RANDOMISING INITIAL VORTEX POSITIONS ----------------------------------------------------------- ##

r = R * sqrt.(rand(NumVortices)); # generates NumVortices random numbers between 0 and 1 that will initialise the radii, takes the square root, and then multiplies by 0.8*R to get the random radial distances for the point vortices
angles = 2*pi*rand(NumVortices); # generates random angles between 0 and 2*pi for the angular locations of the vortices

# vortexPosReal = hcat(r.*cos.(angles), r.*sin.(angles)) # generates random positions of vortices in 2D Cartesian coordinates from the random radii and angles generated above

A = zeros(NumVortices, 2*NumTime); # writing initial positions to file

A[:,1:2] = [vortexPosReal[:,1] vortexPosReal[:,2]];



## ----------------------------------------------------------- SIMULATION LOOP ----------------------------------------------------------- ##

@time begin # time how long the code takes to run

for i in 1:NumTime-1
    global vortexPosReal = vortexPosReal + RK(vortexPosReal, gamma, DampingCoeff, Omega, dt); # does the integration step from the RK function defined above
    
    A[:,(2*i+1):(2*i+2)] = [vortexPosReal[:,1] vortexPosReal[:,2]]; # writes the result to the file
    
end

output_name = raw".\Text file outputs\\" * string(NumVortices) * raw"_vortex_damping_rotation.txt";

writedlm(output_name, A, " "); # writes the output to file

end