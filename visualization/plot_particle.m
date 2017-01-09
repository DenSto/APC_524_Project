% plot particle trajectories
% plots track_<mpi_rank>_<part_rank>.dat output files
% these files contain test particle position and velocity data

function [t,x,v] = plot_particle(mpi_rank,part_rank,do_plot)

if nargin < 3
    do_plot = 1;
end

dname = './';
fname = [dname 'track_' num2str(mpi_rank) '_' num2str(part_rank) '.dat'];

dat=dlmread(fname);

t=dat(:,1);
xx = dat(:,2);
xy = dat(:,3);
xz = dat(:,4);
vx = dat(:,5);
vy = dat(:,6);
vz = dat(:,7);

x=[xx xy xz];
v=[vx vy vz];


if do_plot
    
        % parameters for plots of components
        figure;
        FS = 14;
        LW = 2;
        nrows=2;
        ncols=1;
        
        % plot 3 components of position
        subplot(nrows,ncols,1);
        plot(t,xx,'linewidth',LW); hold all ;
        plot(t,xy,'linewidth',LW);
        plot(t,xz,'linewidth',LW);
        set(gca,'fontsize',FS); box on;
        %     xlabel('Time');
        ylabel('Position');
        legend('x_x','x_y','x_z');
        
        % plot 3 components of velocity
        subplot(nrows,ncols,2);
        plot(t,vx,'linewidth',LW); hold all ;
        plot(t,vy,'linewidth',LW);
        plot(t,vz,'linewidth',LW);
        set(gca,'fontsize',FS); box on ;
        xlabel('Time');
        ylabel('Velocity');
        legend('v_x','v_y','v_z');
        
        suptitle('Components of Position and Velocity');

        % parameters for 3D plots
        FS2 = 10;
        mstyle='.';
        lstyle='-';
        
        % plot position in 3-space
        figure;
        h=plot3(xx,xy,xz);
        set(h,'marker',mstyle,'linestyle',lstyle,'color','k');
        set(gca,'fontsize',FS2); box on;
        xlabel('x'); ylabel('y'); zlabel('z');
        title('Particle Trajectory in Space');
    
        % plot velocity in 3-space
        figure;
        h=plot3(vx,vy,vz);
        set(h,'marker',mstyle,'linestyle',lstyle,'color','k');
        set(gca,'fontsize',FS2); box on;
        xlabel('v_x'); ylabel('v_y'); zlabel('v_z');
        title('Particle Trajectory in Phase Space');
        
end


end