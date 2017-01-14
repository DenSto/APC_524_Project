% gets indices of unique values in array vals
% helper function for plot_field
%
% input:
% vals is an array
%
% output:
% dex is an array of indices such that vals(dex) does not have
% any consecutive values equal
%
% creates files:
% none

function dex = unique_dex(vals)

n = numel(vals);

dup_dex = [];
for i = 1:n-1
    if vals(i) == vals(i+1)
        dup_dex = [dup_dex i];
    end
end
dex = setdiff(1:n,dup_dex);

return 

end