using LinearAlgebra
using DelimitedFiles
using BenchmarkTools

## --------------------------------------------------------------- FUNCTIONS FOR RK IMPLEMENTATION --------------------------------------------------------------- ##

## this function calculates r_{ij}^2 from an array of complex numbers that store the coordinates of each point vortex ##
function r_squared(pos_Real, i, j)

    return (real(pos_Real[i]) - imag(pos_Real[j]))^2 + (real(pos_Real[i]) - imag(pos_Real[j]))^2

end


## this function calculates and returns the complex number that stores the x- and y-coordinates of the ith image vortex - the real component of the output gives the x-coordinate while the imaginary part gives the y-coordinate ##
function posImage(pos_Real, i, R)

    return (R^2)*pos_Real[i]/abs2(pos_Real[i])

end


## this function returns and calculates the right-hand side of the ODE for x and y (with damping) in an array with each row representing a different point vortex and the first and second columns representing the x- and y-coordinates, respectively ##
function derivative(pos_Real, gamma, damping_Coeff, Omega)

    n = size(pos_Real,1)
    vortexVelocity = zeros(ComplexF64,n,1)

    z_ij = 1 ./ (pos_Real' .- pos_Real); # calculates the 2D array that will store the values of every 1/(z_i - z_j) - note that the diagonal entries will be NaN
    z_ij[1:n+1:n^2] .= 0; # sets the diagonal entries equal to zero so that we can do a matrix mulitiplication instead of the summation in the equations of motion

    # this does the summation in the equations of motion but vectorised
    vortexVelocity = (1im/2*pi)*z_ij'*gamma;

    # this part adds the term that accounts for the rotating frame to the local superfluid velocity
    vortexVelocity = vortexVelocity - 1im*Omega*conj(pos_Real)

    #= # this part will first calculate the local superfluid velocity, v_i, for each point vortex location, and then will add the correction for damping
    for i in 1:n
        ## THE COMMENTED OUT STUFF IN THIS LOOP IS FOR THE BOUNDARY CONDITIONS, WHICH WE WILL IGNORE FOR THE SIMULATION
        for j in 1:n
            if j==i
                continue
               # vec[i,1] = vec[i,1] + (1/2*pi)*gamma[j]*(pos_Real[i,2] - posImage(pos_Real,j,R)[2])/r_bar_squared(pos_Real,i,j) # calculates the x-component
               # vec[i,2] = vec[i,2] - (1/2*pi)*gamma[j]*(pos_Real[i,1] - posImage(pos_Real,j,R)[1])/r_bar_squared(pos_Real,i,j) # calculates the y-component
            else
                vortexVelocity[i] = vortexVelocity[i] + (1im/2*pi)*gamma/(pos_Real[j] - pos_Real[i]) # + (1/2*pi)*gamma[j]*(pos_Real[i,2] - posImage(pos_Real,j,R)[2])/r_bar_squared(pos_Real,i,j) # calculates the x-component
            end
        end
        
    end =#

    #= # this part adds the correction from the rotating frame to the local superfluid velocity
    superfluidVelocity[:,1] = superfluidVelocity[:,1] - Omega*pos_Real[:,2]
    superfluidVelocity[:,2] = superfluidVelocity[:,2] + Omega*pos_Real[:,1] =#

    # this part calculates the equations of motion for the actual vortices from the local superfluid velocity calculated above
    vortexVelocity = (1 - 1im*damping_Coeff)*vortexVelocity;

    return conj(vortexVelocity) # returns the velocity of each vortex as an array of complex numbers, each element of the array corresponding to a different vortex - the real part corresponds to the x-coordinates and the imaginary part corresponds to the y-coordinate

end


## this function does the integration step for Runge-Kutta - calculates the part we need to add to each poisition to calculate the positions at the next time step (this will return an array, with each row for a different point vortex, and then the two columns represent x and y) ##
function RK(pos_Real, gamma, damping_Coeff, Omega, dt)
    
    # calculates f0 as just the right-hand side of the ODE
    f0 = derivative(pos_Real, gamma, damping_Coeff, Omega)

    # calculates f1, which is the right-hand side of the ODE at the current position plus 0.5*dt*f0
    f1 = derivative(pos_Real + 0.5*dt*f0, gamma, damping_Coeff, Omega)

    # calculates f2, which is the right-hand side of the ODE at the current position plus 0.5*dt*f1
    f2 = derivative(pos_Real + 0.5*dt*f1, gamma, damping_Coeff, Omega)

    # calculates f3, which is the right-hand side of the ODE at the current position plus 0.5*dt*f1
    f3 = derivative(pos_Real + dt*f2, gamma, damping_Coeff, Omega)

    return (dt/6)*(f0 + 2*f1 + 2*f2 + f3)

end


## --------------------------------------------------------------- DEFINING CONSTANTS AND INITIALISING VARIABLES --------------------------------------------------------------- ##

NumVortices = 1000; # this is the number of vortices in the simulation
dt = 0.001; # time step
tmax = 20; # maximum integration time
NumTime = floor(Int,tmax/dt); # this is the number of time steps

R = 1.0; # radius of superfluid

gamma = ones(ComplexF64, (NumVortices, 1)); # circulation of each point vortex is generated as being equal to 


DampingCoeff = 0.3; # damping coefficient
Omega = 5; # angular velocity of the rotating frame

#global vortexPosReal = zeros(ComplexF64,NumVortices,1); # initialises array that will be used to store the initial (complex) point vortex positions

## ----------------------------------------------------------- RANDOMISING INITIAL VORTEX POSITIONS ----------------------------------------------------------- ##

r = 0.8 * R * sqrt.(rand(ComplexF64, NumVortices)); # generates NumVortices random numbers between 0 and 1 that will initialise the radii, takes the square root, and then multiplies by 0.8*R to get the random radial distances for the point vortices
angles = 2*pi*rand(ComplexF64, NumVortices); # generates random angles between 0 and 2*pi for the angular locations of the vortices

global vortexPosReal = r.*cos.(angles) + r.*sin.(angles)im; # generates random positions of vortices in 2D Cartesian coordinates from the random radii and angles generated above

A = zeros(NumVortices, 2*NumTime); # writing initial positions to file

A[:,1:2] = [real(vortexPosReal[:]) imag(vortexPosReal[:])];


## ----------------------------------------------------------- SIMULATION LOOP ----------------------------------------------------------- ##

@time begin # time how long the code takes to run

for i in 1:NumTime-1
    global vortexPosReal = vortexPosReal + RK(vortexPosReal, gamma, DampingCoeff, Omega, dt);

    A[:,(2*i+1):(2*i+2)] = [real(vortexPosReal[:]) imag(vortexPosReal[:])];
    

end

output_name = raw".\Text file outputs\\" * string(NumVortices) * raw"_vortex_damping_rotation.txt";

writedlm(output_name, A, " ") # writes the output to file

end