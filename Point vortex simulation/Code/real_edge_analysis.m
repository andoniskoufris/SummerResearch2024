% this function will take in the array with all the positions in complex, coordinates and then return R vs theta for each vortex on the edge
function array = real_edge_analysis(positions, numEdgeVortices, radius)
    
    Nt = size(positions, 1); % this extracts the number of time steps from the number of rows in the positions array

    array = zeros(2*Nt, numEdgeVortices); % initialises the array that will return theta against R for the vortices on the edge - these come in pairs of rows, the first being the angle and the second being the radius
    
    array(1:2:2*Nt - 1, :) = angle(positions(:, 1:numEdgeVortices)); % fills in the first row of each pair of rows with the angles of the vortices on the edge - the vortices on the edge are in the first columns of the position array
    
    array(2:2:2*Nt, :) = abs(positions(:, 1:numEdgeVortices)) - radius; % fills in the second row of each pair of rows with the radii of the vortices on the edge (minus the radius of the stationary state) - the vortices on the edge are in the first columns of the position array
    
    % this will sort the data in each pair of rows in increasing phase so that the plot looks correct and doesn't connect back up on itself
    for j=1:Nt
        
        [B, index] = sort(array(2*j-1,:));
        
        array(2*j-1,:) = array(2*j-1,index);
        array(2*j,:) = array(2*j,index);

    end

end