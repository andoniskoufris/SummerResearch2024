% this function will take in the positions array after I have used ode45 with the single vortex and returns the Fourier coefficients with the corresponding wave numbers
function array = fourier_edge_analysis(positions, numEdgeVortices, radius)
    
    Nt = size(positions, 1); % this extracts the number of time steps from the number of rows in the positions array

    edge_real_coefficients = real_edge_analysis(positions, numEdgeVortices, radius); % gets the angles and radii of each vortex in pairs of rows
    
    k = (0:(numEdgeVortices-1))/numEdgeVortices; % defines the default array of k wave numbers that the nufft function will use
    
    array = zeros(2*Nt, numEdgeVortices); % initialises the array that will return the fft of theta against R for the vortices on the edge - these come in pairs of rows, the first being the angle and the second being the radius
    
    % fills in the first row of each pair with the wave numbers that the nufft function will use
    for j=1:2:Nt-1
        array(j, :) = k(:); 
    end
    
    % fills in the second row of each pair of rows with the Fourier transform of the radial data - I need to do this in a for loop because it doesn't like when I use colon stuff for nufft
    for j=1:Nt
        array(2*j, :) = nufft(edge_real_coefficients(2*j-1, :), edge_real_coefficients(2*j, :));
    end

end