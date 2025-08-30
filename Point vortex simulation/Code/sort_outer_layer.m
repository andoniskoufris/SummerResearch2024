function positions_sorted = sort_outer_layer(positions, radius, tolerance)
    
    edge_vortices = identify_outer_layer(positions, radius, tolerance); % give this function the complex coordinates of the vortices and returns an array with just the vortices in the edge layer, sorted by 
    
    positions_sorted = zeros(1, size(positions,2)); % generates a row vector that will store the positions of the new vortex positions, the first entries of which will be the vortices on the outer layer so that we can easily find those vortices that are were originally in the outer layer and how they change over time
    
    positions_sorted(1, 1:size(edge_vortices, 2)) = edge_vortices(:);  % sets the first n entries of the sorted array equal to the positions of the vortices on the edge of the stationary state, so we can easily keep track of them as we keep track of the excitations in the edge layer

    positions_sorted(1, size(edge_vortices, 2)+1:size(positions_sorted,2)) = positions(size(positions,1), (1-tolerance)*radius >= abs(positions(:))); % sets the rest of the entries in the array equal to all the other vortices not in the edge

end