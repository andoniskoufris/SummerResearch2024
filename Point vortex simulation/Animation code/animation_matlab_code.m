clf


data=importdata("C:\Users\skouf\Documents\University\2024\Research\Simulations\Point vortex simulation\julia_test.txt");

radius = 1.1 * max(data(size(data,1), :));

v = VideoWriter('test_thingy','MPEG-4');
open(v);

x0 = 600;
y0 = 400;

width = 450;
height = 400;

num_particles = size(data, 2)/2;

figure(1)
plot(data(1, 1:2:(num_particles-1)), data(1, 1:2:num_particles),'.','MarkerSize', 14, 'Color','black')
%viscircles([0,0],30,'Color','k')
hold off
axis([-radius radius -radius radius])
set(gcf, "Position", [x0,y0,width,height]);
f = getframe(gcf);

%[im,map] = rgb2ind(f.cdata,256,'nodither');
%im(:,:,1,1)=rgb2ind(f.cdata,map,'nodither');
%imwrite(im,map,'julia_test_animation.gif','gif','LoopCount',Inf,'DelayTime',0)

writeVideo(v,f)


for j=2:round(size(data, 1)) - 0.0*round(size(data, 1))
    
    if mod(j,20) == 0
        clf
        plot(data(1, 1:2:(num_particles-1)), data(1, 1:2:num_particles),'.','MarkerSize', 14, 'Color','black')
        %viscircles([0,0],30,'Color','k')
        axis([-radius radius -radius radius])
        
        set(gcf, "Position", [x0,y0,width,height]);
        f = getframe(gcf);

        writeVideo(v,f)
        pause(0.01)

        % [im,map] = rgb2ind(f.cdata,256,'nodither');
        % f = getframe;
        % imwrite(im,map,'julia_test_animation.gif','gif','WriteMode','append','DelayTime',0.04)

    end
    
end
close(v)
