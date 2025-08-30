clear all; close all

%% ----------------------------------------------------------- DEFINING CONSTANTS ----------------------------------------------------------- %%

dt = 0.01; % time step
tmax = 25; % maximum integration time
NumTime = floor(tmax/dt); % this is the number of time steps

tspan = gpuArray(double(linspace(0,tmax,NumTime)));

R = 5.0; % radius of superfluid


DampingCoeff = 0.9; % damping coefficient
Omega = 20; % angular velocity of the rotating frame


%% ----------------------------------------------------------- RANDOMISING INITIAL VORTEX POSITIONS ----------------------------------------------------------- %%

% ----------------------------------------------- UNIFORM TRIANGULAR LATTICE INITIAL CONDITION ----------------------------------------------- %

vortexPosReal = make_triangular_lattice(0.1, R); % generates the initial condition from the make_triangular_lattice.m function

NumVortices = size(vortexPosReal,1); % gets the number of particles from the array created from the triangular lattice function if this is the initial condition I'm using
gamma = ones(NumVortices, 1); % circulation of each point vortex is generated as being equal to 1

vortexPosReal = gpuArray(double(vortexPosReal));
gamma = gpuArray(double(gamma));


%% ----------------------------------------------------------- SIMULATION LOOP FOR REACHING STATIONARY STATE ----------------------------------------------------------- %%

tic

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t_stationary,z_stationary] = ode45(@(t_stationary,z_stationary) derivative_Matt(t_stationary,z_stationary,gamma,DampingCoeff,Omega), tspan, vortexPosReal, options);


%% ----------------------------------------------------------- SIMULATION LOOP FOR ADDING VORTICES TO BOUNDARY ----------------------------------------------------------- %%

R0 = max([max(real(z_stationary(NumTime,:))),max(imag(z_stationary(NumTime,:))),-min(real(z_stationary(NumTime,:))),-min(imag(z_stationary(NumTime,:)))]); % calculates the approximate radius of the stationary state by finding the maximum of the bounds of the stationary state in the x- and y-directions

tspan_perturb = gpuArray(double(linspace(0, tmax, 0.1*NumTime))); % array of time steps for the simulation after we have introduced the blob to perturb the system at the boundary

NumTime_perturb = size(tspan_perturb,2); % gets the number of time steps from the time span over which we simulate the perturbed system


%% ----------------------------------------------------------- TRACKS THE EDGE LAYER VORTICES ----------------------------------------------------------- %%

Tolerance = 0.018;

edge_vortices = identify_outer_layer(z_stationary(NumTime,:), R0, Tolerance); % finds the vortices in the outer edge layer by finding those whose radii are within a tolerance defined above of the distance in the second argument I've given to the function

NumEdgeVortices = size(edge_vortices, 1); % finds the number of vortices on the edge from the size of the array above that finds the vortices in the edge layer

z_sorted = sort_outer_layer(z_stationary(NumTime,:), R0, Tolerance); % sorts the vortices in the stationary, with the first of the vortices in the array being those vortices on the edge


% -------------------------------------------------------- CREATE A UNIFORMLY RANDOM CIRCULAR BLOB TO PERTURB THE BOUNDARY -------------------------------------------------------- %

% -------------------------------------------------------- CREATE PERTURBING VORTICES WITH POSITION BY HAND -------------------------------------------------------- %

NumPerturbationVortices = 1; % this is the number of vortices we will add to "perturb" the boundary

vorticesBlob = zeros(NumPerturbationVortices,1); % initialises array that will be used to store the perturbing vortex positions

angle_edge_vortices = angle(edge_vortices); % finds the phase of each edge vortex - the identify_outer_layer function sorts the positions by phase

[theta_0, index_0] = min(abs(angle(edge_vortices) - 0)); % finds the angle and index of the angle of the edge vortices closest to zero

theta_1 = (angle_edge_vortices(index_0 + 1) + angle_edge_vortices(index_0)) / 2; % angle that the perturbing vortex will make with the positive x-axis
% theta_2 = angle_edge_vortices(index_0 + 1);

vorticesBlob = 1.01 *R0*exp(1i*theta_1);

stationary_perturbation_init = zeros(NumVortices + NumPerturbationVortices, 1); % initialises the (complex) array that will store the initial conditions for after we perturb the boundary of the system in the stationary state

% % fills in the array that will store the initial condition with both the vortices in the stationary state and the perturbing vortices
% stationary_perturbation_init(1:NumVortices) = z_stationary(NumTime, :);
% stationary_perturbation_init(NumVortices+1:NumVortices + NumPerturbationVortices) = vorticesBlob;
% stationary_perturbation_init = gpuArray(double(stationary_perturbation_init));

% does things with the sorted array with the initial conditions
stationary_perturbation_init(1:NumVortices) = z_sorted(:);
stationary_perturbation_init(NumVortices+1:NumVortices + NumPerturbationVortices) = vorticesBlob;
stationary_perturbation_init = gpuArray(double(stationary_perturbation_init));


gamma_perturb = gpuArray(double(ones(NumVortices + NumPerturbationVortices,1))); % creates an array that stores the circulation for each vortex with the additional perturbing vortices

[t_perturb, z_perturb] = ode45(@(t_perturb,z_perturb) derivative_Matt(t_perturb,z_perturb,gamma_perturb,0,Omega), tspan_perturb, stationary_perturbation_init, options); % does the simulation for the system that starts in the stationary state with the added perturbing vortices
% 
% grid_space_coefficients = real_edge_analysis(z_perturb, size(edge_vortices, 1), R0); % gets the real space coefficients of angle versus theta after adding the single perturbation outside the edge of the system
% fourier_space_coefficients = fourier_edge_analysis(z_perturb, size(edge_vortices, 1), R0); % gets the Fourier space coefficients and wave numbers after adding the single perturbation outside the edge of the system
% 
% % writes the real space and the fourier space coefficients both to files
% writematrix(real(grid_space_coefficients), strcat("C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\", num2str(NumVortices + NumPerturbationVortices,'%.0f'), "_real_coefficient.txt"), 'Delimiter', " ") % writes the real space coefficients to file from analysing edge modes
% writematrix(real(fourier_space_coefficients), strcat("C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\", num2str(NumVortices + NumPerturbationVortices,'%.0f'), "_fourier_coefficient.txt"), 'Delimiter', " ") % writes the fourier space coefficients to file from analysing edge modes


% -------------------------------------------------------- FUNCTION THAT FINDS THE OUTER BOUNDARY LAYER AND THEN ALSO PERTURBS IT -------------------------------------------------------- %

% outer_layer_vortices = identify_outer_layer(z_stationary(NumTime,:), R0, 0.018); % defines an array with just the positions of the outer layer vortices I find with my function
% new_edge_vortices = perturb_outer_layer(outer_layer_vortices, 5, R0);
% 
% full_density_perturbation = gpuArray(double( density_perturbation(z_stationary(NumTime,:), R0, 0.018, 5) ));
% 
% %plot(real(z_stationary(NumTime,:)), imag(z_stationary(NumTime,:)), '.k', 'MarkerSize',6)
% hold on
%     %plot(real(outer_layer_vortices), imag(outer_layer_vortices), '.r', 'MarkerSize',6)
%     %plot(real(new_edge_vortices), imag(new_edge_vortices), '.b', 'MarkerSize',6)
%     plot(real(full_density_perturbation), imag(full_density_perturbation), '.b', 'MarkerSize',6)
% hold off
% 
% [t_perturb, z_perturb] = ode45(@(t_perturb,z_perturb) derivative_Matt(t_perturb,z_perturb,gamma,0,Omega), tspan_perturb, full_density_perturbation, options); % does the simulation for the system that starts in the stationary state with the density perturbation


%% ----------------------------------------------------------- WRITE OUTPUT TO FILE ----------------------------------------------------------- %%


% writing perturbation from the addition of vortices to file
A_perturb = zeros(NumVortices + NumPerturbationVortices, 2*NumTime_perturb); % initialises array that will be written to file for the simulation after the boundary has been "perturbed" with the additional vortices

% this loop and following code writes the positions after the system has been perturbed to a file
for j=1:NumTime_perturb
    % writes positions of vortices from the simulation after we have perturbed the boundary to array that will be written to file
    A_perturb(:,2*j-1) = real(z_perturb(j,:));
    A_perturb(:,2*j) = imag(z_perturb(j,:));
end

output_name_perturb = strcat("/home/s4698066/Point vortex code", num2str(NumVortices + NumPerturbationVortices,'%.0f'), "_single_perturbation.txt");

writematrix(A_perturb, output_name_perturb, 'Delimiter', " ") % writes the output to file
% 
toc


% % writing density perturbation to file
% A_density_perturb = zeros(NumVortices, 2*NumTime_perturb); % initialises array that will be written to file for the simulation after the boundary has been "perturbed" with the additional vortices
% 
% % this loop and following code writes the positions after the system has been perturbed to a file
% for j=1:NumTime_perturb
%     % writes positions of vortices from the simulation after we have perturbed the boundary to array that will be written to file
%     A_density_perturb(:,2*j-1) = real(z_perturb(j,:));
%     A_density_perturb(:,2*j) = imag(z_perturb(j,:));
% end
% 
% output_name_density_perturb = strcat("/home/s4698066/Point vortex code", num2str(NumVortices,'%.0f'), "_density_perturbation.txt");
% 
% writematrix(A_density_perturb, output_name_density_perturb, 'Delimiter', " ") % writes the output to file
% toc


%% ----------------------------------------------------------- FUNCTION DEFINITIONS ----------------------------------------------------------- %%

% this function calculates the velocity of each vortex in complex coordinates and returns it as a vector
function vortexVelocity = derivative(t, pos_Real, gamma, damping_Coeff, Omega)
    
    N = size(pos_Real,1); % extracts the number of vortices from the 
    z_ij = 1./(pos_Real - pos_Real.'); % defines a 2D array that stores the values of 1/(z_i - z_j)
    z_ij(1:N+1:N^2) = 0; % sets the diagonal elements of the array equal to zero - he is indexing by only one number which in MATLAB will go along the row first and then down a column once you've reached the end of the row
    
    vortexVelocity = -1i*Omega*conj(pos_Real) + (1i/2*pi)*z_ij*gamma; % calculates the summation in x- and y-coordinates using complex coordinates - this has been vectorised so that instead of actually doing the summation it's vectorised
    
    %vortexVelocity = vortexVelocity - 1i*damping_Coeff*vortexVelocity; % adds dissipation effect and returns the velocity in complex coordinates - the real part is the velocity in the x-direction and the imaginary part is the velocity in the y-direction

    %vortexVelocity = conj(vortexVelocity); % returns the velocity in complex coordinates - the real part is the velocity in the x-direction and the imaginary part is the velocity in the y-direction
    
    vortexVelocity = conj(vortexVelocity) + 1i*damping_Coeff*conj(vortexVelocity); % adds dissipation effect and returns the velocity in complex coordinates - the real part is the velocity in the x-direction and the imaginary part is the velocity in the y-direction

end


% this function calculates the velocity of each vortex in complex coordinates using the procedure in Matt's code
function vortexVelocity = derivative_Matt(t, pos_Real, gamma, damping_Coeff, Omega)
    
    N = size(pos_Real,1); % extracts the number of vortices from the 
    z_ij = 1./(pos_Real - pos_Real.'); % defines a 2D array that stores the values of 1/(z_i - z_j)
    z_ij(1:N+1:N^2) = 0; % sets the diagonal elements of the array equal to zero - he is indexing by only one number which in MATLAB will go along the row first and then down a column once you've reached the end of the row
    
    vortexVelocity = -conj( 1i*Omega*pos_Real ) - (1i)*z_ij*gamma; % calculates the summation in x- and y-coordinates using complex coordinates - this has been vectorised so that instead of actually doing the summation it's vectorised
    
    vortexVelocity = conj(vortexVelocity) - 1i*damping_Coeff*conj(vortexVelocity); % adds dissipation effect

end