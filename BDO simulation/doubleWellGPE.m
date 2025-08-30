N = 256;
L = 40;
dx = L/N;
dk = 2*pi/L;

x = (-N/2:N/2-1)*dx;
k = (-N/2:N/2-1)*dk;
k = fftshift(k);


V0 = 2.0;
a = 1.5;
V = -V0*(exp(-(x-a).^2) + exp(-(x+a).^2));
g = 2;

psi = exp(-(x-a).^2/2);
psi = psi/sum(abs(psi).^2)/dx;

dt = 1e-2;
Nt = 10000;

for jj =1:Nt
    
    psi = exp(-1i*(V+g*(abs(psi).^2))*dt/2).*psi;
    psi = fft(psi);
    psi = exp(-1i*k.^2*dt/2).*psi;
    psi = ifft(psi);
    psi = exp(-1i*(V+g*(abs(psi).^2))*dt/2).*psi;
    
    if (floor(jj/20) == jj/20)
    figure(1)
    clf
    plot(x,V)
    hold on
    plot(x,abs(psi).^2)
    ylim([-1.5*V0, 1])
    xlim([-1/2 1/2]*L)
    drawnow
    end
    
end
