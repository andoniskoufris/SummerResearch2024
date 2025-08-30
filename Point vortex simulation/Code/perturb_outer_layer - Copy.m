% this function will take in the complex positions of the vortices on the edge of the vortex matter and will return the complex positions of the vortices just now with a slightly different density profile of my choosing
function new_outer_positions = perturb_outer_layer(outer_positions, n, R0)

    phases = angle(outer_positions); % this generates the phase of each vortex on the outer layer whose positions I give in the "outer_positions" array I pass to the function - I add pi just so that the angles are between 0 and 2*pi rather than -pi and pi
    radii = abs(outer_positions); % this generates the phase of each vortex on the outer layer whose positions I give in the "outer_positions" array I pass to the function
    
    index_random = randi([1 + n, size(outer_positions,2) - n]); % chooses a random vortex number to centre our perturbation around

    vortex_phase_centre = phases(index_random); % chooses the phase of a random vortex on the edge to clump the distribution

    phase_min = phases(index_random - n); % finds the minimum angle that we will distribute the random points between that will make up the disturbance

    phase_max = phases(index_random + n); % finds the maximum angle that we will distribute the random points between that will make up the disturbance
    
    phase_range = abs(phase_max - phase_min)/2; % finds the range of phases we will be distributing the new vortices that follow our distribution of choice over
    
    angle_range = phase_range; % converts the range of angles for the vortices we will be replacing into a physical length along the perimeter of the stationary state by just multiplying the angular interval by the radius of the stationary state

    % random_points = rand(n,1); % generates the random points we feed into the distribution functions below
    random_points = rand(1,n);

    % --------------------------------------------- GENERATES RANDOM ANGLES ACCORDING TO DISTRIBUTION OF CHOICE --------------------------------------------- %

    A = 1/pi;
    
    new_phases = zeros(1, 2*n); % initialises array that will store the new vortex positions on either side of the central point

    % Lorentzian distribution attempting to use inverse sampling method
    new_phases(1:n) = vortex_phase_centre + angle_range*rejection_sampling(n,A); % generates the new phases of the vortex positions that will replace the old ones about some central point
    new_phases(n+1:2*n) = vortex_phase_centre - angle_range*rejection_sampling(n,A); % generates the new phases of the vortex positions that will replace the old ones about some central point
    
    % % attempted Lorentzian distribution by hand
    % hand_angles = (1/0.66) * [0.02; 0.08; 0.17; 0.32; 0.6]; % hand chosen angles for a Lorentzian profile
    % 
    % new_phases(1:n) = vortex_phase_centre - angle_range*hand_angles; % generates the new phases of the vortex positions that will replace the old ones about some central point
    % new_phases(n+1:2*n) = vortex_phase_centre + angle_range*hand_angles; % generates the new phases of the vortex positions that will replace the old ones about some central point
    
    
    % ------------------------------ DENSITY PERTURBATION ------------------------------ %
    % replacement_positions(1:n) = radii(index_random - n:index_random - 1).*exp(1i*new_phases(1:n)); % generates the 2*n+1 vortex positions that will replace the current vortex positions
    % 
    % replacement_positions(n+1:2*n) = radii(index_random + 1:index_random + n).*exp(1i*new_phases(n+1:2*n)); % generates the 2*n+1 vortex positions that will replace the current vortex positions as being alll at the same radius R0 and at the angles generated above
    % 
    % new_outer_positions = outer_positions; % sets the array that will contain all the new vortex positions equal to the array with the old vortex positions
    % 
    % new_outer_positions(index_random-n:index_random-1) = replacement_positions(1:n); 
    % new_outer_positions(index_random+1:index_random+n) = replacement_positions(n+1:2*n);
    
    % ------------------------------ RADIUS PERTURBATION ------------------------------ %
    
    % perturbation_height = 0.03*R0; % sets the height of the perturbation
    % 
    % radius_perturbation = R0 + perturbation_height*A./((exp(2)*(phases(index_random - n:index_random + n) - vortex_phase_centre)/angle_range).^2 + A^2); % perturbs the radius of the edge vortices according to a Lorentzian distribution
    % 
    % new_outer_positions = outer_positions; % sets the array that will contain all the new vortex positions equal to the array with the old vortex positions
    % 
    % new_outer_positions(index_random - n:index_random + n) = radius_perturbation.*exp(1i*phases(index_random - n:index_random + n));

    % % ------------------------------ RADIUS AND DENSITY PERTURBATION ------------------------------ %

    % perturbation_height = 0.03*R0; % sets the height of the perturbation
    % 
    % replacement_positions = zeros(2*n,1);
    % 
    % radius_perturbation = R0 + perturbation_height*A./((exp(2)*(new_phases - vortex_phase_centre)/angle_range).^2 + A^2); % perturbs the radius of the edge vortices according to a Lorentzian function
    % 
    % replacement_positions(1:n) = radius_perturbation(1:n).*exp(1i*new_phases(1:n)); % generates the 2*n+1 vortex positions that will replace the current vortex positions
    % replacement_positions(n+1:2*n) = radius_perturbation(n+1:2*n).*exp(1i*new_phases(n+1:2*n));
    % 
    % new_outer_positions = outer_positions; % sets the array that will contain all the new vortex positions equal to the array with the old vortex positions
    % 
    % new_outer_positions(index_random-n:index_random-1) = replacement_positions(1:n); % replaces the positions of the vortices currently on the edge with those whose density has been perturbed, as well as the radii
    % new_outer_positions(index_random+1:index_random+n) = replacement_positions(n+1:2*n);
    % new_outer_positions(index_random) = 1.03*new_outer_positions(index_random); % moves the middle vortex in the perturbation out by the appropriate amount
    
    % ------------------------------ DENSITY PERTURBATION ON TOP OF EXISTING VORTICES ------------------------------ %
    
    new_positions(1:n) = radii(index_random - n:index_random - 1).*exp(1i*new_phases(1:n)); % generates the 2*n+1 vortex positions that will replace the current vortex positions

    new_positions(n+1:2*n) = radii(index_random + 1:index_random + n).*exp(1i*new_phases(n+1:2*n)); % generates the 2*n+1 vortex positions that will replace the current vortex positions as being alll at the same radius R0 and at the angles generated above

    new_outer_positions = outer_positions; % sets the array that will contain all the new vortex positions equal to the array with the old vortex positions

    new_outer_positions = [new_outer_positions, new_positions]; % appends the new vortex positions to the current vortex positions

end