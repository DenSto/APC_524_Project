% plot a single plane of the field data at specified times
%
% input:
% field is a [nx,ny,nz,nt] size array of field data as a function
% of space and time, such as created by plot_field.m
%
% slicedir is an integer determining which direction to slice the field in
% slicedir = 1 is constant x, 2 is constant y, 3 is constant z
% default: slicedir = 3
% 
% islice is an integer index to slice, between 1 and size(field,slicedir)
% default: islice = 1
%
% itimes is an integer array specifying which time slices to plot
% set itimes to a negative scalar to plot all available time slices
% default: itimes = -1
% 
% do_save = 1 saves plots to file
% do_save = 0 omits saving
% default: do_save = 1
%
% returns:
% none
%
% creates files (if do_save == 1):
% if numel(itimes) > 1 or itimes > 0
% field.png 2D plot of field 
% if numel(itimes) > 0 or itimes < 0
% field_mov.avi 2D movie of field vs time

function [] = slice_field(field,slicedir,islice,itimes,do_save)

% set default values
if nargin < 5
    do_save = 1;
end

if nargin < 4
    itimes = -1;
end

if nargin <3
    islice = 1;
end

if nargin < 2
    slicedir = 3;
end

% get field dimensions
nt=size(field,4);
nx=size(field,1);
ny=size(field,2);
nz=size(field,3);

% set times for movie
if numel(itimes) ==1 && itimes < 0
    itimes=1:nt;
end
nt = numel(itimes);

% set to 1 to make a contour plot
% set to 2 to make a surface plot
plot_type = 2;

if nt > 1
    do_mov = 1;
    % for saving movie
    F(nt) = struct('cdata',[],'colormap',[]);
    fsave='field_mov.avi';
else
    do_mov = 0;
    fsave='field.png';
end

f = figure;
set(f,'color','w'); 
for i = itimes
    if (slicedir == 1)
        dat = reshape(field(islice,:,:,i),[ny nz]);
        datmin = min(min(min(field(islice,:,:,:))));
        datmax = max(max(max(field(islice,:,:,:))));
    elseif (slicedir == 2)
        dat = reshape(field(:,islice,:,i),[nx nz]);
        datmin = min(min(min(field(:,islice,:,:))));
        datmax = max(max(max(field(:,islice,:,:))));
    elseif (slicedir == 3)
        dat = reshape(field(:,:,islice,i),[nx ny]);
        datmin = min(min(min(field(:,:,islice,:))));
        datmax = max(max(max(field(:,:,islice,:))));
    end
    
    if plot_type == 1
        contourf(dat);
    elseif plot_type == 2
        surf(dat);
        zlim([datmin datmax]);
    end
    box on ; colorbar;
    
    drawnow;
    pause(.1);
    
    % save a frame for the video
    if do_mov && do_save
        lighting phong
        set(f,'Renderer','zbuffer');
        F(i) = getframe(f);
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
        save_and_close([fstr '.png'],f,do_save,0)
    end
end

end