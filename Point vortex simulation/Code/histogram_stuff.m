clear all; close all

v = VideoWriter('histogram_edge','MPEG-4');

v.FrameRate = 10;

open(v);

x0 = 500;
y0 = 400;

width = 700;
height = 500;


z_real = readmatrix(strcat("C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\Text file outputs\10981_density_perturbation.txt"));

NumTime = size(z_real,2)/2;
NumVortices = size(z_real,1);

z_complex = zeros(NumTime, NumVortices); % initialises array that will store the complex positions of vortices - each row is a time step, each vortex along the columns

for jj=1:NumTime
    z_complex(jj,:) = z_real(:,2*jj-1) + 1i*z_real(:,2*jj);
end

radius_max = max(abs(z_complex(1,:)));

edge_indices = abs(z_complex(1,:))>0.991*radius_max; % finds indices of vortices in the edge in the stationary state

% plot(real(z_complex(1,:)), imag(z_complex(1,:)), '.')
% hold on
% plot(real(z_complex(1,abs(z_complex(1,:))>0.991*radius_max)),imag(z_complex(1,abs(z_complex(1,:))>0.991*radius_max)),'.')
% hold off

for jj=1:NumTime

    if (floor(jj/20) == jj/20)
    %z_edge = z_complex(jj, abs(z_complex(1,:))>0.991*radius_max); % finds the number of vortices in the outer layer
   

    clf

    hold on
    
    histogram(angle(z_complex(jj,edge_indices)),200);
    %histogram(angle(z_edge),200)

    drawnow

    set(gcf, "Position", [x0,y0,width,height]);
    f = getframe(gcf);
    writeVideo(v,f)
    pause(0.1)
    end
end

close(v)