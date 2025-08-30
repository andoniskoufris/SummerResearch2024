% this function finds the indices of the vortices that satisfy the condition of being within some tolerance from the radius I give the function
function array = identify_outer_layer(positions, radius, tolerance)
    
    array1 = positions((1-tolerance)*radius < abs(positions(:))); % extracts the position of the vortices that satisfy the condition of being within a tolerance (which you give to the function) of the estimated radius of the stationary cluster that I give the function 

    [B, index] = sort(angle(array1) + pi); % obtains the order of the phase of each vortex on the outer layer

    array = array1(index); % sorts the positions of the vortices on the boundary by the phase of the complex coordinates

end