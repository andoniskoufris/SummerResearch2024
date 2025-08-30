N = 256;
L = 20;
dx = L/N;
dk = 2*pi/L;

x = (-N/2:N/2-1)*dx;
k = (-N/2:N/2-1)*dk;
k = fftshift(k);

% defining constants
dbar = 1/(8*pi);
Ubar = 1;
gamma = 10;
A=10;


v = VideoWriter('KDV_equation_animation','MPEG-4');
open(v);

x0 = 500;
y0 = 400;

width = 600;
height = 400;


psi = 8*sech(2*(x+8)).^2 + 2*sech(x+1).^2; % initial condition for the number of vortices along the linearised boundary - this initial condition will be one soliton, as on page 4 of the (pre-print) paper

dt = 0.4/N^2;
Nt = 100000;

figure(1)
plot(x,abs(psi).^2)
ylim([-1.1*max(abs(psi).^2), 1.1*max(abs(psi).^2)])
xlim([-1/2 1/2]*L)

set(gcf, "Position", [x0,y0,width,height]);
f = getframe(gcf);

writeVideo(v,f)


for jj =1:Nt
   
    psi = fft(psi); % computes the FFT of psi that we will use to perform the time-evolution of the analytic part of the KDV
    
    psi = exp(1i*k.^3*dt).*psi; % computes the time-evolution of the analytical term in the KDV
    
    % psi = ifft(psi); % inverse Fourier-transforms our solution so that we can perform Euler's method for the time-stepping of the non-analytical part of the KDV as per the split-step method
    
    psi = ifft( psi - 6*1i*dt*k.*fft(real(ifft(psi)).^2) );
    

    if (floor(jj/100) == jj/100)
    figure(1)
    clf
    hold on
    plot(x,abs(psi).^2)
    ylim([-1.1*max(abs(psi).^2), 1.1*max(abs(psi).^2)])
    xlim([-1/2 1/2]*L)
    drawnow

    set(gcf, "Position", [x0,y0,width,height]);
    f = getframe(gcf);

    writeVideo(v,f)
    pause(0.01)

    end
    
end

close(v)