close all; clear all

tic

set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

N = 4000;
L = 2000;
%L = 30; % for gaussian initial condition
dx = L/N;
dk = 2*pi/L;

x = (-N/2:N/2-1)*dx;
k = (-N/2:N/2-1)*dk;
k = fftshift(k);

k_sign = sign(k); % this provides an array that gives the sign of the wave numbers that we can use in our formula for the Hilbert transform of a Fourier transform

% defining constants
dbar = 1/(8*pi);
Gamma = 1;
Omega = 20;
l = sqrt(Gamma/2*Omega);
Ubar = Gamma/sqrt(16*pi*0.0001);
A1 = 16; % the free parameter for the 1-soliton solution from the paper, which controls the speed, width and height of the wave
A2 = 4; % the free parameter for the second soliton
A3 = 60; % the free parameter for the third soliton

v = VideoWriter('three-soliton-2','MPEG-4');
open(v);

x0 = 500;
y0 = 400;

width = 700;
height = 500;

options = odeset('RelTol',1e-8,'AbsTol',1e-8); % options for using ode45 in the 

num1 = 8*dbar*A1./((x + 750).^2 + A1^2) + 8*dbar*A2./((x + 450).^2 + A2^2) + 8*dbar*A3./((x - 300).^2 + A3^2); % initial condition for the number of vortices along the linearised boundary - this initial condition will be one soliton, as on page 4 of the (pre-print) paper
% num1 = (1/pi)*16./((x+700).^2 + 16^2); % Gaussian initial condition
% num2 = 0.02*exp(-(x+700).^2/(2*30^2));


%dt = 0.4/N^2;
dt = 1e-3;
Nt = 700000;

tspan = linspace(0, dt*Nt, Nt);



figure(1)
%subplot(2,1,1)
plot(x,real(num1))
%xlabel('$x$ (arb. units)')
ylabel('$n$ (dimensionless)')
ylim([-1.1*abs(min(real(num1))), 1.1*abs(max(real(num1)))])
xlim([-1/2 1/2]*L)
%title('Lorentzian initial condition')

% subplot(2,1,2)
% plot(x,real(num2))
xlabel('$x$ (arb. units)')
% ylabel('$n$ (dimensionless)')
% ylim([-1.1*abs(min(real(num2))), 1.1*abs(max(real(num2)))])
% xlim([-1/2 1/2]*L)
% title('Gaussian initial condition')

set(gcf, "Position", [x0,y0,width,height]);
f = getframe(gcf);

writeVideo(v,f)

Ubar = 16;
Gamma = 400;

for jj=1:Nt-1

    % ------------------------------- TIME STEP WITH CUSTOM RK4 FUNCTION ------------------------------- %

    num1 = fft(num1); % computes the FFT of n that we will use to perform the time-evolution of the analytic part of the BDO
    %num2 = fft(num2); % computes the FFT of n that we will use to perform the time-evolution of the analytic part of the BDO

    num1 = exp(-1i*k*Ubar*dt).*num1; % computes the time-evolution of the analytical term in the BDO
    %num2 = exp(-1i*k*Ubar*dt).*num2; % computes the time-evolution of the analytical term in the BDO

    num1 = ifft( num1 + RK4(num1, k, Gamma, dbar, dt) );
    %num2 = ifft( num2 + RK4(num2, k, Gamma, dbar, dt) );

    % ------------------------------- TIME STEP WITH ODE45 ------------------------------- %

    % [t, f] = ode45(@(t,f) derivative(t,f,k,gamma,dbar), tspan(jj:jj+1), vortexPosReal, options); % calculates the derivative 
    % num = ifft( num + f ); % does the full time step with the result from ode45 for a single time step

    if (floor(jj*dt) == jj*dt)
    figure(1)
    clf
    hold on
    %subplot(2,1,1)
    plot(x,real(num1))
    ylim([-1.1*abs(min(real(num1))), 1.1*abs(max(real(num1)))])
    xlim([-1/2 1/2]*L)
    xlabel('$x$ (arb. units)')
    ylabel('$n$ (dimensionless)')
    %title('Lorentzian initial condition')
    
    %subplot(2,1,2)
    % hold on
    % plot(x,real(num2))
    % ylim([-1.1*abs(min(real(num2))), 1.1*abs(max(real(num2)))])
    % xlim([-1/2 1/2]*L)
    % xlabel('$x$ (arb. units)')
    % ylabel('$n$ (dimensionless)')
    % title('Gaussian initial condition')

    drawnow
    
    set(gcf, "Position", [x0,y0,width,height]);

    f = getframe(gcf);

    writeVideo(v,f)
    pause(0.01)

    end

end

close(v)

toc


% This function is my custom-defined Hilbert transform from the definition on Wikipedia, but using the identity stated on Wikipedia for the Hilbert transform of a Fourier transform - this function returns the FFT of the Hilbert transform of a function, so the array that this function returns will need to be inverse FFTed to be plotted
function array = HilbertTransform_identity(c_n, k)

    array = -1i*sign(k).*c_n; % the Hilbert transform of a Fourier transform is just i*sign(k_n)*c_n - this line is for if I just want to give the function k, and then compute the sign_k inside the function each time step instead of just giving it the array that gives the sign of the wave numbers, since this array is never going to change
    
    % array = -1i*sign_k.*c_n; % the same line as before, but for when I  want to speed up my code and just compute sign(k) outside the loop at the start of my code so it doesn't need to be computed at each time step
   
end


% this function calculates and returns the non-analytical part of the BDO equation - takes in n in Fourier space and returns it in Fourier space
function array = nonAnalyticalBit(f, k, gamma, dbar)

    array = 0.5*gamma*1i*k.*fft(ifft(f).^2) + 2*gamma*dbar*k.^2.*HilbertTransform_identity(f,k); 

end


% this function does the RK4 time step with a custom-defined function
function array = RK4(f, k, gamma, dbar, dt)

    f0 = nonAnalyticalBit(f, k, gamma, dbar);

    f1 = nonAnalyticalBit(f + 0.5*dt*f0, k, gamma, dbar);

    f2 = nonAnalyticalBit(f + 0.5*dt*f1, k, gamma, dbar);

    f3 = nonAnalyticalBit(f + dt*f2, k, gamma, dbar);

    array = (dt/6)*(f0 + 2*f1 + 2*f2 + f3);

end