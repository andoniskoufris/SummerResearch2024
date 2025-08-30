function perturbed_positions = density_perturbation(positions, radius, tolerance, n)
    
    % --------------------------------------------- THIS IS STUFF FROM THE IDENTIFY_OUTER_LAYER FUNCTION --------------------------------------------- %    
    
    array1 = positions((1-tolerance)*radius < abs(positions(:))); % extracts the position of the vortices that satisfy the condition of being within a tolerance (which you give to the function) of the estimated radius of the stationary cluster that I give the function 

    [B, index] = sort(angle(array1) + pi); % obtains the order of the phase of each vortex on the outer layer

    array = array1(index); % sorts the positions of the vortices on the boundary by the phase of the complex coordinates

    % --------------------------------------------- THEN PERTURB OUTER LAYER USING FUNCTION --------------------------------------------- %
    
    perturbed_boundary_layer = perturb_outer_layer(array, n, radius); % perturbs the outer layer using the perturb_outer_layer function

    % --------------------------------------------- THEN PUT THESE POSITIONS BACK INTO FULL POSITION ARRAY --------------------------------------------- %
    
    perturbed_positions = zeros(1, size(positions,2)+2*n); % generates a row vector that will store the positions of the new vortex positions after we have perturbed the density of the boundary 
    
    perturbed_positions(1:size(array, 2)+2*n) = perturbed_boundary_layer(:); 

    perturbed_positions(size(array, 2)+2*n+1:size(perturbed_positions,2)) = positions((1-tolerance)*radius >= abs(positions(:)));

end