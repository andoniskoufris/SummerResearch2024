clear
clc
%close all


Nvec = [2 5 10 20 50 100 200 500 1000 2000]
%for kk =1:length(Nvec)



for kk = 1%:1000
N = 331; % number of particles (I think)
R = 1; % radius of initial configuration
Omega = N/R^2; % angular velocity of the rotating superfluid
gamma0 = 0.1; % 
rho0 = N/pi/R^2;
Rt = @(tt) R*sqrt(1 + 2*pi*rho0*gamma0*tt);
kappa = ones(1,N);
theta = linspace(0,2*pi,100);

%for run = 1
%disp(num2str(run))
z0 = 0.5*R*sqrt(rand(N,1)).*exp(1i*2*pi*rand(N,1));
z0 = makeHexagonal();
z0 = z0/max(abs(z0))*0.9*exp(1i*0.01)
%z0 = R*(randn(N,1) + 1i*randn(N,1))/sqrt(2);
M = sum(abs(z0).^2) /N;
z0 = z0./sqrt(M)*sqrt(0.5);
%z0 = 0.25*(randn(N,1) + 1i*randn(N,1));

figure(1)
clf
plot(z0,'.')
xlim([-1 1]*R*1.2)
ylim(xlim)
hold on
plot(exp(1i*theta),'--k')

%%
Nt = 10000;
tf = 10/Omega/gamma0;
ss = @(tt) log(1+rho0*gamma0*2*pi*tt)/R^2/rho0/gamma0/2/pi;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
s = linspace(0,tf,Nt);
%s = ss(t);
tic
[t,z] = ode45(@(t,z) PointVortexPlane(t,z,kappa,gamma0,N,rho0,R),s,z0,options);
toc
z_save(:,kk) = z(end,:);
u_save(:,kk) = getVelocity(z_save(:,kk),kappa,N);
%%

zbin = 0:0.005:1;
for jj = 1:length(zbin)-1
ind = z_save(:) > zbin(jj) & z_save(:) < zbin(jj+1);
U(jj) = mean(abs(u_save(ind)));
U2(jj) = mean(angle(u_save(ind)));
end

% rad(kk) = max(abs(z(end,:)))
% E(kk) = Energy(z(end,:),N)/N^2
% M(kk) = mean(abs(z(end,:)).^2)
figure(3)
plot(z(end,:),'.r')

drawnow
disp(num2str(kk))


figure(4)
clf
[counts,vals] = hist(abs(z_save(:)),200);
counts = counts./vals/sum(counts)/(vals(2)-vals(1));
plot(vals,counts)

yyaxis right
plot(zbin(1:end-1),(U))
hold on
plot(zbin(1:end-1),U2)
drawnow
end
% h = figure(3);
% axis square
% for jj = 1:1:Nt
%     cla
%     %plot(z(jj,:)/Rt(t(jj))*(1+rho0*gamma0*2*pi*t(jj))^(-1i*Omega/gamma0/2/pi/rho0),'.')
%     plot(z(jj,:),'.r')
%     hold on
%     plot(exp(1i*theta),'-k')
%     xlim([-1 1]*1.5)
%     ylim(xlim)
%     title(num2str(jj))
%     drawnow
%     F(jj) = getframe(gcf);
%     E(jj) = Energy(z(jj,:),N);
% end

% figure(10)
% hold on
% plot(Omega*t,E/N^2)

%%
% % create the video writer with 1 fps
%   writerObj = VideoWriter('vortexmatter.avi');
%   writerObj.FrameRate = 60;
%   % set the seconds per image
%     % open the video writer
%     open(writerObj);
%     % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

%end

function Z = makeHexagonal()
base=1;
X=[];
Y=[];
for num=1:10
    x=zeros(num*6,1);
    y=zeros(num*6,1);
    x(1:6)=base*num*cos(2*pi/6.*(0:5));
    y(1:6)=base*num*sin(2*pi/6.*(0:5));
    if num>1
        for q=1:num-1
           start_x=x(2)-q*base;
           radi0=sqrt(start_x^2+y(2)^2);
           start_alpha=1/3*pi+pi*1/3*1/(num)*q;
           x(q*6+1:q*6+6)=radi0*cos(start_alpha+pi/3.*(1:6));
           y(q*6+1:q*6+6)=radi0*sin(start_alpha+pi/3.*(1:6));
        end
    end
    X=[X; x; ];
    Y=[Y; y; ];
    
end
X = [X;  0];
Y = [Y;  0];
scatter(X,Y)

Z = X+1i*Y;
end 
    