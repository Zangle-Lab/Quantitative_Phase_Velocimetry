function [ip, jp,results] = findvalley_v3_rev4(A)
%function [ip, jp] = findpeak_v2(A)
%function to find the peak of A, using subpixel grid Xsub,Ysub, of size
%gs x gs
%v2 uses gaussian fitting to attempt subpixel accuracy
% sz = size(A);
% gs = 3; 
% [Xsub,Ysub] = meshgrid(1:gs,1:gs); %using 3 x 3 grid to find answer with subpixel accuracy
%v3 uses fast gaussian ps function from matlab file exchange

% A = uint16(A.*1000); %speed up by casting to int? didn't work in practice..
%%
%old way, also requires pad array
% [~, km] = max(A(:));
% [im,jm] = ind2sub(sz,km);

% B = padarray(A,[(gs-1)./2,(gs-1)./2]); %pad array so we dont run into problems at the edge
% Bsub = B(im-(gs-1)/2+1:im+(gs-1)/2+1,jm-(gs-1)/2+1:jm+(gs-1)/2+1);

% %use weighted centroids:
% x0 = (sum(Xsub(:).*Bsub(:))./sum(Bsub(:)));
% y0 = (sum(Ysub(:).*Bsub(:))./sum(Bsub(:)));
%%

sz = size(A);
[row,col] = find(ismember((A), min(min(A(:)))));
if length(row)==1 && length(col)==1
    if row>1 && row<sz(1)-1
        if col>1 && col<sz(2)-1
            E=A(row-1:row+1,col-1:col+1);
            a=max(E(:));
            Q=double(a-E);
            results = psfFit(Q);
%             ip = (((3-results(1))/5)-(((sz(1)+1)/2)-col));
%             jp = (((3-results(2))/5)-(((sz(2)+1)/2)-row));
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
% if length(row)==1 && length(col)==1
%     dely=(W(row+1,col)-W(row-1,col));
%     delx=(W(row,col+1)-W(row,col-1));
%     ip=col-((sz(1)+4+1)/2)+delx;
%     jp=row-((sz(1)+4+1)/2)+dely;
% else
%     ip=0;
%     jp=0;
% end
% if length(row)==1 && length(col)==1
%     E=fliplr(flipud(W(row-2:row+2,col-2:col+2)));
%     results = psfFit(E);
%     ip = -(((3-results(1))/5)+(((sz(1)+5)/2)-col));
%     jp = -(((3-results(2))/5)+(((sz(2)+5)/2)-row));
% else
%     ip=0;
%     jp=0;
% end
% opts.positive = true;
% results = autoGaussianSurf(Xsub,Ysub,A,opts);
% ip = -results.x0; %amount to shift in x so that next lines up with current
% jp = -results.y0; %amount to shift in y so that next lines up with current
