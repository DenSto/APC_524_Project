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