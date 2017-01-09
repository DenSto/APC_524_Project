% plot particle trajectories

function [Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,rho] = plot_fields(fname,do_plot)

if nargin < 2
    do_plot = 1;
end

if nargin < 1
    dname='./';
    fname=[dname 'output.h5'];
end

datpath='/fields/';

% read each field from the data file
Ex = h5read(fname,[datpath 'Ex']);
Ey = h5read(fname,[datpath 'Ey']);
Ez = h5read(fname,[datpath 'Ez']);

Bx = h5read(fname,[datpath 'Bx']);
By = h5read(fname,[datpath 'By']);
Bz = h5read(fname,[datpath 'Bz']);

Jx = h5read(fname,[datpath 'Jx']);
Jy = h5read(fname,[datpath 'Jy']);
Jz = h5read(fname,[datpath 'Jz']);

rho = h5read(fname,[datpath 'rho']);

% later: read the physical grid from the data file 

if do_plot
    FS = 14;
    
    figure;
    [x,y,z] = meshgrid(-2:.2:2,-2:.25:2,-2:.16:2);
    v = x.*exp(-x.^2-y.^2-z.^2);
    xslice = [-1.2,.8,2];
    yslice = 2;
    zslice = [-2,0];
    slice(x,y,z,v,xslice,yslice,zslice)
    set(gca,'fontsize',FS); box on; 
    
end


end