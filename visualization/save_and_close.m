% helper function to close or save figures 
%
% input: 
% fname is a string specifying file name to save to 
% do_save = 1 to save the figure to file, 0 to omit saving 
% do_close = 1 to close the figure after saving, 0 to leave open 
% default behavior is do_save = 1, do_close = 1
% 
% returns: 
% none 


function save_and_close(fname,fig,do_save,do_close)
if nargin < 3 
    do_close = 1; 
end 

if nargin < 2
    do_close = 1; 
end 

if do_save
    hgexport(fig, fname, hgexport('factorystyle'), 'Format', 'png');
end
if do_close
    close(fig);
end
end