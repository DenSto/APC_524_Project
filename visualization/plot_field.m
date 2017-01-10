% plot particle trajectories

% function [Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,rho] = plot_fields(fname,do_plot)
function [t,x,y,z,field] = plot_field(fstr,do_plot,do_movie)

% plot by default
if nargin < 2
    do_plot = 1;
end

% expect output file output.h5 to be in this directory
% if elsewhere, change dname to path
dname='./';
fname=[dname 'output.h5'];

% path within hdf5 file
datpath='/fields/';

% set which field to read
datset=[datpath fstr];
try
    field = h5read(fname,datset);
catch
    error('Error: h5read error, probably specified string incorrectly');
end

x = h5readatt(fname,datset,'x');
y = h5readatt(fname,datset,'y');
z = h5readatt(fname,datset,'z');

nx=numel(x);
ny=numel(y);
nz=numel(z);

% set the time to an index for now
% later time will be written as well
nt = size(field,1);
t = 1:nt;

% for some reason, fields are written in t,z,y,x ordering
% permute into more natural x,y,z,t
field = permute(field,[4 3 2 1]);

xmin=min(x); xmax=max(x);
ymin=min(y); ymax=max(y);
zmin=min(z); zmax=max(z);

% verify that x,y,z grid matches field size
if (nx ~= size(field,1))
    error('error: x grid does not match size of field');
elseif (ny ~= size(field,2))
    error('error: y grid does not match size of field');
elseif (nz ~= size(field,3))
    error('error: z grid does not match size of field');
end


% order must be reversed because of meshgrid quirks
% not possible to use ndgrid with slice function
[X,Y,Z] = meshgrid(y,x,z);

if do_plot
    
    FS = 14;
    
    % set to 1 to plot the faces of the box
    % set to 0 to turn them off
    do_faces=0;
    
    % set to 1 to plot slices through the middle of the domain
    % set to 0 to turn them off
    % (do_mid and do_faces both = 1 will show faces)
    do_mid=1;
    
    % set to 1 to make a symmetric colorscale around 0
    % set to 0 to use min and max values for colorscale
    % (0 means the middle value in the colorscale corresponds to
    % the average value, not zero)
    norm_color = 0;
    
    xslice = [];
    yslice = [];
    zslice = [];
    if do_faces
        xslice = [xslice xmin xmax];
        yslice = [yslice ymin ymax];
        zslice = [zslice zmin zmax];
    end
    
    if do_mid
        xmean = mean([xmin xmax]);
        ymean = mean([ymin ymax]);
        zmean = mean([zmin zmax]);
        
        xslice = [xslice xmean];
        yslice = [yslice ymean];
        zslice = [zslice zmean];
    end
    
    % set the colormap for plotting
    % other options include cmap = 'parula','jet',etc
    cmap = jet(1024);
    
    figure;
    
    if do_movie
        times = 1:nt;
    else
        times = floor(nt/2); % arbitrary time
    end
    
    % loop over all times
    for it = times
        v = reshape(field(:,:,:,it),[nx ny nz]);
        
        % set bounds for colorscale
        minv = min(min(min(v)));
        maxv = max(max(max(v)));
        
        if minv == maxv
            eps = 1e-10;
            minv = minv - eps;
            maxv = maxv + eps;
        end
        
        if norm_color
            cbound = max([abs(minv) maxv]);
            cvec = cbound*[-1 1];
        else
            cvec = [minv maxv];
        end
        
        % plot 3D data slice
        % (note: similarly here, order of xslice and yslice are
        % counterintuitive but I think this is the correct ordering)
        slice(X,Y,Z,v,yslice,xslice,zslice)
        set(gca,'fontsize',FS); box on;
        colormap(cmap); colorbar;
        try
            caxis(cvec);
        catch
            disp(num2str(cvec));
        end
        
        pause(.05);
    end
    
end


end