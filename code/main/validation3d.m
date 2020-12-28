clc
clear all
close all

% Reading the results from FORTRAN code
fileID1=fopen('animation3d_parallel.dat','r');

dum=fscanf(fileID1,'%s',1);
npoints=fscanf(fileID1,'%d',1);
dum=fscanf(fileID1,'%s',1);
ntime=fscanf(fileID1,'%d',1);

fprintf("Npoints:%d, ntime:%d\n",npoints,ntime)

Nx=sqrt(npoints);
Ny=Nx;

tt=cputime;

Nfile=Nx*Ny; %Number of lines in the input file
dt=1e-6;
Nt=ntime;

dum=fscanf(fileID1,['%f %f'],[2 Nx*Ny]);
xFOR=dum(1,:)';
yFOR=dum(2,:)';

dx=0.5/(Nx-1);
dy=0.5/(Ny-1);

loc=zeros(Nx*Ny,2);
for iy=1:Ny
    for ix=1:Nx
        loc((iy-1)*Nx+ix,1)=round(xFOR((iy-1)*Nx+ix)/dx)+(Nx+1)/2;
        loc((iy-1)*Nx+ix,2)=round(yFOR((iy-1)*Nx+ix)/dy)+(Ny+1)/2;
    end
end

xFOR=-0.25:dx:0.25;
yFOR=-0.25:dy:0.25;
[XX YY]=meshgrid(xFOR,yFOR);

F(ntime) = struct('cdata',[],'colormap',[]);
for it=1:ntime
    time=fscanf(fileID1,['%f'],[1]);
    dum=fscanf(fileID1,['%f'],[1 Nx*Ny]);
    
     
    disp=zeros(Nx,Ny);
    for iy=1:Ny
    for ix=1:Nx
        indx=loc((iy-1)*Nx+ix,1);
        indy=loc((iy-1)*Nx+ix,2);
        disp(indx,indy)=dum((iy-1)*Nx+ix);
    end
    end
   
   contourf(XX,YY,disp)
   colorbar
%  plot(xFOR,disp(:,50))
% ylim([-5e-9 5e-9])
title(['Parallel t=',num2str(time)])

    drawnow
    F(it) = getframe;
end

%%

% for it=1:200
%     it
%     time=fscanf(fileID1,['%f'],[1]);
%     dum=fscanf(fileID1,['%f'],[1 Nx*Ny]);
%     
%     if(mod(it,20)==0)
%         
%     disp=zeros(Nx,Ny);
%     for iy=1:Ny
%     for ix=1:Nx
%         indx=loc((iy-1)*Nx+ix,1);
%         indy=loc((iy-1)*Nx+ix,2);
%         disp(indx,indy)=dum((iy-1)*Nx+ix);
%     end
%     end
%    
%     figure
%    contourf(XX,YY,disp)
%    colorbar
% title(['Parallel t=',num2str(time)])
% 
%      end
% end
% 
% 
