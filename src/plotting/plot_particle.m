% plot particle trajectories

function [t,x,v] = plot_particle(mpi_rank,part_rank,do_plot)

if nargin < 3
    do_plot = 1;
end

fname = ['track_' num2str(mpi_rank) '_' num2str(part_rank) '.dat'];
f = fopen(fname);

nt = 1e5;
t = zeros(nt,1);
x = zeros(nt,3);
v = x;

fmt='%f';
iter = 0;
while ~feof(f)
    iter = iter + 1;
    
    tl = fscanf(f,fmt,1);
    xl = fscanf(f,fmt,3);
    vl = fscanf(f,fmt,3);
    
    if ~isempty(tl)
        t(iter) = tl;
        x(iter,:) = xl;
        v(iter,:) = vl;
    else
        break;
    end
end

fclose(f);

if do_plot
    FS = 14;
    LW = 2;
    figure;
    plot(t,x(:,1),'linewidth',LW); hold on ;
    plot(t,x(:,2),'linewidth',LW);
    plot(t,x(:,3),'linewidth',LW);
    set(gca,'fontsize',FS);
    xlabel('Time');
    ylabel('Position');
    legend('x','y','z');
    box on;
    
    figure;
    plot(t,v(:,1),'linewidth',LW); hold on ;
    plot(t,v(:,2),'linewidth',LW);
    plot(t,v(:,3),'linewidth',LW);
    set(gca,'fontsize',FS);
    xlabel('Time');
    ylabel('Velocity');
    legend('x','y','z');
    box on;
    
end


end