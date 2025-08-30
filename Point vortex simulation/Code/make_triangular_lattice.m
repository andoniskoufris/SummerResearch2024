function array = make_triangular_lattice(h_dist,R)

% Triangular grid information - h_dist defines horizontal spacing between points in triangular lattice
v_dist = sqrt(h_dist^2-(h_dist/2)^2); % defines vertical spacing between points in triangular lattice

% This is the extent of our triangular lattice in each direction
x_lim = 2*R;
y_lim = x_lim;

% this code and the following while loop generates the points in the grid - this grid will be square but the spacing between points will be triangular
trigrid = [];
y_current = 0;
xx = 0;
displacement = 0;

while y_current < y_lim
    if displacement == 0
        xx = [0:h_dist:x_lim]';
        yy = ones(length(xx), 1)*y_current;
        displacement = 1;
    else
        xx = [h_dist/2:h_dist:x_lim]';
        yy = ones(length(xx), 1)*y_current;
        displacement = 0;
    end
    trigrid = [trigrid; [xx,yy]];
    y_current = y_current + v_dist;
end

trigrid = trigrid - (x_lim/2);

% this following loop will go through the points from the triangular lattice made above in the shape of a square, and will get rid of points that are outside of a circle of a certain radius from the centre

j=1;

while j<size(trigrid,1)
    if sqrt(trigrid(j,1)^2 + trigrid(j,2)^2)>x_lim/2
        trigrid(j,:) = [];
    else
        j = j+1;
    end
end

trigrid(size(trigrid,1),:) = [];

array = zeros(size(trigrid,1),1);

for j=1:size(trigrid,1)
    array(j) = trigrid(j,1) + 1i*trigrid(j,2);
end

end

% figure()
% plot(trigrid(:,1), trigrid(:,2), 'o', 'markersize', 2);
% grid on;
% xlim([-x_lim-h_dist, x_lim+h_dist]);
% ylim([-y_lim-v_dist, y_lim+v_dist]);
% 
% writematrix(trigrid,'triangular_lattice.txt','Delimiter',' ');