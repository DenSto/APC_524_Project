% plot single field in 3D at specified times
% plots output.h5 output files
% output.h5 files contain field data
% assumes output.h5 is in current working directory
%
% input:
% fstr is a string specifying which field to plot. Options are:
% Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho
% though output.h5 will only contain the fields specified in input.txt
%
% do_plot = 1 to make plots,
% do_plot = 0 to omit plotting (useful for reading hdf5 into matlab variables)
% default: do_plot = 1
%
% itimes is an integer array specifying which time slices to plot
% set itimes to a negative scalar to plot all available time slices
% default: itimes = -1
%
% do_save = 1 saves plots to file
% do_save  0 omits saving
% default: do_save = 1
%
% returns:
% t is a vector of times of length nt when field was recorded
% x is a vector with the x grid points of the field, same for y,z
% field is an array of size (nx,ny,nz,nt) with the full field
% data at each point in space and time
%
% creates files (if do_save == 1):
% if numel(itimes) > 1 or itimes > 0
% <fstr>.png 3D plot of field
% if numel(itimes) > 0 or itimes < 0
% <fstr>_mov.avi 3D movie of field vs time

function [t,x,y,z,field] = plot_field(fstr,do_plot,itimes,do_save)

% plot and save by default
if nargin < 4
    do_save = 1;
end

if nargin < 3
    itimes = -1;
end

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

% for some reason, fields are written in t,z,y,x ordering
% permute into more natural x,y,z,t
field = permute(field,[4 3 2 1]);

x = h5readatt(fname,datset,'x');
y = h5readatt(fname,datset,'y');
z = h5readatt(fname,datset,'z');

% remove any redundancies in the field data 
% (not usually an issue)
xdex = unique_dex(x); 
ydex = unique_dex(y); 
zdex = unique_dex(z); 

x = x(xdex); 
y = y(ydex);
z = z(zdex); 

field = field(xdex,ydex,zdex,:); 

nx=numel(x);
ny=numel(y);
nz=numel(z);

% set the time to an index for now
% later time will be written as well
t = h5read(fname,'/time');
nt = numel(t);

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
elseif (nt ~= size(field,4))
    error('error: t grid does not match size of field');
end

% order must be reversed because of meshgrid quirks
% not possible to use ndgrid with slice function
[X,Y,Z] = meshgrid(y,x,z);

if do_plot
    
    if numel(itimes)==1 && itimes < 0
        itimes=1:nt;
    end
    
    if numel(itimes) > 1
        do_mov = 1;
    else
        do_mov = 0;
    end
    
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
    norm_color = 1;
    
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
    %     cmap = make_rbcmap(1024);
    %     cmap = jet(1024);
    cmap = 'parula';
    
    f = figure;
    
    if do_mov && do_save
        F(nt)= struct('cdata',[],'colormap',[]);
        fsave=[fstr '_mov.avi'];
    end
    
    % loop over all times
    for it = itimes
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
        
        % choose your own slices to add
        % set physical values, not slice indices
        %         xslice=[xslice 3];
        %         yslice=[yslice 2.53];
        %         zslice=[zslice .93];
        
        % plot 3D data slice
        % (note: similarly here, order of xslice and yslice are
        % counterintuitive but I think this is the correct ordering)
        slice(X,Y,Z,v,yslice,xslice,zslice)
        set(gca,'fontsize',FS); box on; colorbar;
        
        % this should be fine but has some strange consequences
        % for matlab 2013
        colormap(cmap);
        colorbar;
        caxis(cvec);
        
        % labels also dumb due to meshgrid
        xlabel('y'); ylabel('x'); zlabel('z');
        axis tight; axis equal;
        
        drawnow; 
        pause(.05);
        
        % save a frame for the video
        if do_mov && do_save
            lighting phong
            set(f,'Renderer','zbuffer');
            F(it) = getframe(f);
        end
    end
    
end

if do_save
    if do_mov
        v=VideoWriter(fsave,'Motion JPEG AVI');
        v.FrameRate=60;
        open(v);
        writeVideo(v,F);
        close(v);
    else
        save_and_close([fstr '.png'],f,do_save,~do_plot)
    end
end

end