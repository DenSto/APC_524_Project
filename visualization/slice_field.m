% plot a single slice of the field

% slicedir = 1,2,3
function [] = slice_field(field,slicedir,islice,it)
if nargin < 4
    it = -1;
end

if nargin <3
    islice = 2;
end

if nargin < 2
    slicedir = 3;
end

nt=size(field,4);
nx=size(field,1);
ny=size(field,2);
nz=size(field,3);

if it < 0
    it=1:nt;
end

nt = numel(it);

if nt > 1
    % for saving movie
    F(nt) = struct('cdata',[],'colormap',[]);
    fsave='field_mov.avi';
end

f = figure;
for i = it
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
    
    %     contourf(dat);
    surf(dat);
    zlim([datmin datmax]);
    box on ; colorbar;
    
    if nt > 1
        drawnow;
        pause(.1);
        lighting phong
        set(f,'Renderer','zbuffer')
        F(i) = getframe(f);
    end
    
end

if nt > 1
    v=VideoWriter(fsave,'Motion JPEG AVI');
    v.FrameRate=8;
    open(v);
    writeVideo(v,F);
    close(v);
end

end