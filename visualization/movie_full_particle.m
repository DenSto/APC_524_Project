% similar to movie_particle except captures 4 subplots simultaneously: 
% 1. x,y,z components as a function of t (1D plot) 
% 2. position(t) in space (3D plot) 
% 3. x,y,z components of velocity 

function [t,x,v] = movie_full_particle(mpi_rank,part_rank,do_save)
if nargin < 3
    do_save = 1;
end

dname = './';
fname = [dname 'track_' num2str(mpi_rank) '_' num2str(part_rank) '.dat'];

dat=dlmread(fname,'',1,0);

t=dat(:,1);
xx = dat(:,2);
xy = dat(:,3);
xz = dat(:,4);
vx = dat(:,5);
vy = dat(:,6);
vz = dat(:,7);

x=[xx xy xz];
v=[vx vy vz];

nt = numel(t);

xxstr='x'; xystr='y'; xzstr='z';
vxstr='v_x'; vystr='v_y'; vzstr='v_z';
xtstr='Space'; vtstr = 'Phase Space';

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
nparts = 1; 

% size of most recent and least recent marker plotted
msmin = 1;
msmax = 10;

% number of steps to iterate forward
%(e.g. dt = 3 plots on steps 1, 4, 7,etc)
dt=5;

% set to 1 to keep the axis limits fixed for the whole movie
% set to 0 to dynamically adjust them
do_fixed_axes = 0;

%%% end plotting options %%%

% colormap for plotting
movcmap=parula(nhist);
% marker sizes for plotting
msizes=linspace(msmin,msmax,nparts);
% static limits of volume for plotting
eps=.001;
lfrac=1-eps; rfrac=1+eps;
xxl = lfrac*min(xx); xxr = rfrac*max(xx);
xyl = lfrac*min(xy); xyr = rfrac*max(xy);
xzl = lfrac*min(xz); xzr = rfrac*max(xz);
vxl = lfrac*min(vx); vxr = rfrac*max(vx);
vyl = lfrac*min(vy); vyr = rfrac*max(vy);
vzl = lfrac*min(vz); vzr = rfrac*max(vz);

if do_save
    F(nt) = struct('cdata',[],'colormap',[]);
    fsave='full_mov.avi';
end

f=figure('units','normalized','position',[0 .05 1 .85]) ;
set(f,'color','w');

LWt = 1;
LWtraj = 2;
cgrey = .9*ones(1,3);

xt = subplot(2,2,1);
plot_components(t,x);
xt_ylim = ylim;
vt = subplot(2,2,3);
plot_components(t,v);
vt_ylim = ylim;

xtraj = subplot(2,2,2);
vtraj = subplot(2,2,4);

% loop over times
for i=1:dt:nt
    jmin=max([i-nhist+1 1]);
    
    subplot(xt);
    plot_components(t,x,LWt);
    
    subplot(vt)
    plot_components(t,v,LWt);
    
    
    if jmin < i
        subplot(xt);
        tfill = t([jmin i i jmin]);
        xtfill = [xt_ylim(1) xt_ylim(1) xt_ylim(2) xt_ylim(2)];
        fill(tfill,xtfill,cgrey);
        
        subplot(vt);
        vtfill = [vt_ylim(1) vt_ylim(1) vt_ylim(2) vt_ylim(2)];
        fill(tfill,vtfill,cgrey);
    else
        subplot(xt);
        plot(t(i)*[1 1],[xt_ylim(1) xt_ylim(2)],'k');
        subplot(vt);
        plot(t(i)*[1 1],[vt_ylim(1) vt_ylim(2)],'k');
    end
    
    % loop over only the most recent nhist points in color
    for j=jmin:i
        ic = i-j+1;
        jprev = max([1 j-1]);
        
        % plot the previous step with varying color
        subplot(xtraj);
        plot3(xx(jprev:j),xy(jprev:j),xz(jprev:j),'color',movcmap(ic,:),'linewidth',LWtraj);
        % plot nparts most recent markers
        if i-j < nparts
            h=plot3(xx(j),xy(j),xz(j));
            set(h,'marker','o','markerfacecolor',movcmap(ic,:),'markeredgecolor','k','markersize',msizes(nparts-(i-j)));
        end
        hold on ;
        
        subplot(vtraj);
        plot3(vx(jprev:j),vy(jprev:j),vz(jprev:j),'color',movcmap(ic,:),'linewidth',LWtraj);
        % plot nparts most recent markers
        if i-j < nparts
            h=plot3(vx(j),vy(j),vz(j));
            set(h,'marker','o','markerfacecolor',movcmap(ic,:),'markeredgecolor','k','markersize',msizes(nparts-(i-j)));
        end
        hold on ;
    end
    
    % emphasize components actively being plotted 
    subplot(xt); 
    plot_components(t(jmin:i),x(jmin:i,:),LWtraj); 
    axis tight; 
    subplot(vt); 
    plot_components(t(jmin:i),v(jmin:i,:),LWtraj); 
    axis tight; 
    
    % plot full history without color
    if show_full_hist
        jlast = min([jmin nt]);
        subplot(xtraj);
        plot3(xx(1:jlast),xy(1:jlast),xz(1:jlast),'k','linewidth',1);
        subplot(vtraj);
        plot3(vx(1:jlast),vy(1:jlast),vz(1:jlast),'k','linewidth',1);
    end
    
    % plot limits
    if do_fixed_axes
        subplot(xtraj);
        xlim([xxl xxr]); ylim([xyl xyr]); zlim([xzl xzr]);
        subplot(vtraj);
        xlim([vxl vxr]); ylim([vyl vyr]); zlim([vzl vzr]);
    end
    
    % plot labels
    subplot(xt);
    xlabel('Time'); ylabel('Position'); title('Position');
    legend('x','y','z'); 
    box on ; hold off;
    
    subplot(vt);
    xlabel('Time'); ylabel('Velocity'); title('Velocity');
    legend('v_x','v_y','v_z'); 
    box on ; hold off;
    
    subplot(xtraj);
    xlabel(xxstr); ylabel(xystr); zlabel(xzstr);
    title(['Trajectory in ' xtstr ' at t = ' num2str(t(i),'%.2f')]);
    box on ; hold off;
    
    subplot(vtraj);
    xlabel(vxstr); ylabel(vystr); zlabel(vzstr);
    title(['Trajectory in ' vtstr ' at t = ' num2str(t(i),'%.2f')]);
    box on ; hold off;
    
    pause(.001);
    drawnow;
    
    % save a frame for the video
    if do_save
        lighting phong
        set(f,'Renderer','zbuffer');
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