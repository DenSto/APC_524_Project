% plot particle trajectories

function [Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,rho] = read_fields(fname)

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

end