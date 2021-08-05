function [ip, jp,results] = findvalley_v3_rev4(A)

sz = size(A);
[row,col] = find(ismember((A), min(min(A(:)))));
% lowest point in 3 by 3 pixel region around lowest pixel
if length(row)==1 && length(col)==1
    if row>1 && row<sz(1)-1
        if col>1 && col<sz(2)-1
            E=A(row-1:row+1,col-1:col+1);
            a=max(E(:));
            Q=double(a-E);
            results = psfFit(Q);        % Gaussian fitting
            ip=(-((sz(1)+1)/2)+col)-(2-results(1));
            jp=(-((sz(2)+1)/2)+row)-(2-results(2));
        else
            ip=0;
            jp=0;
        end
    else
        ip=0;
        jp=0;
    end
else
    ip=0;
    jp=0;
end

