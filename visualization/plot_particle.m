% plot particle trajectories
% plots track_<mpi_rank>_<part_rank>.dat output files
% track.dat files contain test particle position and velocity data
% assumes track.dat file is in current working directory
%
% input:
% do_plot = 1 to make plots, do_plot = 0 to omit plotting
% (do_plot = 0 simply returns values read from the file)
% do_save = 1 to save plots as .png, do_save = 0 to omit saving
% note: do_save takes precedence over do_plot, such that if do_save = 1
% and do_plot = 0, figures will be created and then immediately closed
% after saving.
%
% returns:
% t is a vector of times of length nt when particle trajectory was recorded
% x is a matrix of size [nt 3] with x,y,z components of position
% v is a matrix of size [nt 3] with x,y,z components of velocity

function [t,x,v] = plot_particle(mpi_rank,part_rank,do_plot,do_save)

if nargin < 4
    do_save = 1;
end
if nargin < 3
    do_plot = 0;
end

dname = './';
fname = [dname 'track_' num2str(mpi_rank) '_' num2str(part_rank) '.dat'];

% second line to deal with new header in track.dat files
% dat=dlmread(fname);
dat = dlmread(fname,'',1,0); 

t=dat(:,1);
xx = dat(:,2);
xy = dat(:,3);
xz = dat(:,4);
vx = dat(:,5);
vy = dat(:,6);
vz = dat(:,7);

x=[xx xy xz];
v=[vx vy vz];

if (do_save || do_plot)
    do_close = ~do_plot;
     
    % parameters for plots of components
    fcomp = figure;
    FS = 14;
    nrows=2;
    ncols=1;
    
    % set to 1 to subtract mean from components (easier to compare) 
    % set to 0 to plot raw values
    do_norm = 0;
    
    if do_norm 
        for i=1:3
            x(:,i) = x(:,i)-mean(x(:,i)); 
            v(:,i) = v(:,i)-mean(v(:,i)); 
        end
    end
    
    % plot 3 components of position
    subplot(nrows,ncols,1);
    plot_components(t,x);
    set(gca,'fontsize',FS); box on;
    %     xlabel('Time');
    ylabel('Position');
    legend('x_x','x_y','x_z');
    
    % plot 3 components of velocity
    subplot(nrows,ncols,2);
    plot_components(t,v);
    set(gca,'fontsize',FS); box on ;
    xlabel('Time');
    ylabel('Velocity');
    legend('v_x','v_y','v_z');
    
    suptitle('Components of Position and Velocity');
    save_and_close('xv_comps.png',fcomp,do_save,do_close);
    
    % parameters for 3D plots
    FS2 = 10;
    mstyle='.';
    lstyle='-';
    
    % plot position in 3-space
    fx = figure;
    h=plot3(xx,xy,xz);
    set(h,'marker',mstyle,'linestyle',lstyle,'color','k');
    set(gca,'fontsize',FS2); box on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Particle Trajectory in Space');
    save_and_close('x_traj.png',fx,do_save,do_close);
    
    % plot velocity in 3-space
    fv = figure;
    h=plot3(vx,vy,vz);
    set(h,'marker',mstyle,'linestyle',lstyle,'color','k');
    set(gca,'fontsize',FS2); box on;
    xlabel('v_x'); ylabel('v_y'); zlabel('v_z');
    title('Particle Trajectory in Phase Space');
    save_and_close('v_traj.png',fv,do_save,do_close)
end

end