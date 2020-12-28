clc
clear all
close all

dt=1e-6;

fileID=fopen('animation2d_parallel.dat' ,'r');

dum=fscanf(fileID,'%s',1);
npoints=fscanf(fileID,'%d',1);
dum=fscanf(fileID,'%s',1);
Tpoints=fscanf(fileID,'%d',1);

xFOR=fscanf(fileID,'%f',npoints);
figure 

F(Tpoints) = struct('cdata',[],'colormap',[]);
for i=1:Tpoints

    time=fscanf(fileID,'%f',1);
    disp=fscanf(fileID,'%f',npoints);

    xlabel('x')
    plot(xFOR,disp)

    title(['t=',num2str(time)])
    ylim([-1.5e-8 1.5e-8])
    drawnow
    F(i) = getframe;

    
end
