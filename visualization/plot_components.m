% helper function called by particle plotting routines
% including: plot_particle.m, movie_particle.m, movie_full_particle.m
% makes 1D plots of components of vectors
%
% input: 
% t is a vector of length nt of times to plot 
% x is a nt x 3 matrix with columns corresponding to x,y,z components 
% of some vector such as position or velocity 
% LW is an optional parameter to set the line width of the lines
% default value is LW = 2
% 
% returns: 
% none 
%
% creates files: 
% none

function [] = plot_components(t,x,LW)

if nargin < 3
    LW = 2;
end

xx = x(:,1);
xy = x(:,2);
xz = x(:,3);

plot(t,xx,'r','linewidth',LW); hold all ;
plot(t,xy,'g','linewidth',LW);
plot(t,xz,'b','linewidth',LW);
end