% makes a movie of particle trajectory in space or phase space
%
% input:
% t is a vector of times of length nt
%
% vec is matrix of size [nt 3] representing x,y,z coordinates of some
% time evolving vector. Anticipated to be position (for plot in real space)
% or velocity (for plot in velocity phase space), but does not have to be.
%
% movtype = 1 for plots of position, movtype = 2 for plots of velocity
% (movtype only affects plot labels and saved file names)
%
% do_save = 1 to save a .avi movie file, do_save = 0 to not save
% default: do_save = 1
%
% returns:
% none
% 
% creates files (if do_save == 1):
% if movtype == 1: 
% x_mov.avi 3D movie of particle position trajectory in space vs time 
% if movtype == 2: 
% v_mov.avi 3D movie of particle velocity trajectory in phase space vs time

function [] = movie_particle(t,vec,movtype,do_save)

% save by default
if nargin < 3
    do_save = 1;
end

% total number of data points
nt = numel(t);

%%% plotting options %%%

% set to 1 to plot black line for all past points
% set to 0 to only plot the nhist most recent points (set below)
show_full_hist=1;

% number of recent points of trajectory to show in color
% set to 0 to plot no points in color (faster)
% nhist = floor(nt/100);
nhist=10;

% number of recent points of trajectory to show with markers
% nparts = floor(nhist/10);
nparts = 0; 

% size of most recent and least recent marker plotted
msmin = 1;
msmax = 10;

% number of steps to iterate forward
%(e.g. dt = 3 plots on steps 1, 4, 7,etc)
dt=2;

% set to 1 to keep the axis limits fixed for the whole movie
% set to 0 to dynamically adjust them
do_fixed_axes = 0;

%%% end plotting options %%%

% labels corresponding to plot of position or velocity
if movtype == 1 % meat space
    tstr = 'Space at t = ';
    xstr='x';
    ystr='y';
    zstr='z';
    savstr='x';
elseif movtype == 2 % phase space
    tstr = 'Phase Space at t = ';
    xstr = 'v_x';
    ystr = 'v_y';
    zstr = 'v_z';
    savstr='v';
end

% coordinates to plot
xx=vec(:,1);
xy=vec(:,2);
xz=vec(:,3);

% colormap for plotting
movcmap=parula(nhist);
% marker sizes for plotting
msizes=linspace(msmin,msmax,nparts);
% static limits of volume for plotting
eps=.001;
lfrac=1-eps; rfrac=1+eps;
xl = lfrac*min(xx); xr = rfrac*max(xx);
yl = lfrac*min(xy); yr = rfrac*max(xy);
zl = lfrac*min(xz); zr = rfrac*max(xz);

if do_save
    F(nt) = struct('cdata',[],'colormap',[]);
    fsave=[savstr '_mov.avi'];
end

f=figure;
set(f,'color','w'); 
% loop over all times
for i=1:dt:nt
    jmin=max([i-nhist+1 1]);
    % loop over only the most recent nhist points in color
    for j=jmin:i
        ic = i-j+1;
        jprev = max([1 j-1]);
        % plot the previous step with varying color
        plot3(xx(jprev:j),xy(jprev:j),xz(jprev:j),'color',movcmap(ic,:),'linewidth',2);
        
        % plot nparts most recent markers
        if i-j < nparts
            h=plot3(xx(j),xy(j),xz(j));
            set(h,'marker','o','markerfacecolor',movcmap(ic,:),'markeredgecolor','k','markersize',msizes(nparts-(i-j)));
        end
        hold on ;
    end
    % plot full history without color
    if show_full_hist
        jlast = min([jmin nt]);
        plot3(xx(1:jlast),xy(1:jlast),xz(1:jlast),'k');
    end
    hold off
    
    % plot labels and limits
    if do_fixed_axes
        xlim([xl xr]); ylim([yl yr]); zlim([zl zr]);
    end
    xlabel(xstr); ylabel(ystr); zlabel(zstr);
    title(['Trajectory in ' tstr num2str(t(i),'%.2f')]);
    box on;
    
    pause(.001);
    drawnow;
    
    % save a frame for the video
    if do_save
        lighting phong
        set(f,'Renderer','zbuffer')
        F(i) = getframe(f);
    end
end

% construct the movie file
if do_save
    F=F(1:dt:nt);
    v=VideoWriter(fsave,'Motion JPEG AVI');
    v.FrameRate=60;
    open(v);
    writeVideo(v,F);
    close(v);
end
end