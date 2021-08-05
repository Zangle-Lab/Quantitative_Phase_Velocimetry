function [SS,BWmap,B1] = imagebackground_poly4_kmeans(DD)
%% kmeans method to segment the cells from background
sz=size(DD);
c = kmeans(DD(:), 3, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
for tt=1:3
   Waq=(CC==tt);
   waq(tt)=mean(mean(Waq.*DD));
end
[row,col] = find(ismember((waq), min(min(waq(:)))));
BWmap1=CC~=col;
seD = strel('diamond',2);
BWfinal = imerode(BWmap1,seD);
BWfinal = bwareaopen(BWfinal, 1000); %remove regions smaller than 10 pixels
se90 = strel('line', 28, 90);
se0 = strel('line', 28, 0);
BWmap=imdilate(BWfinal,[se0 se90]);

%% polyfit using higher order polynomial
IList = DD(~BWmap);
sz = size(DD);
[x,y] = meshgrid(1:sz(2), 1:sz(1));
XList = x(~BWmap);
YList = y(~BWmap);
if sum(~isnan(BWmap)) ~=0
    CFit = polyfitn([XList,YList], IList, 6);
    B1 =(reshape(polyvaln(CFit, [x(:), y(:)]), sz(1), sz(2)));
else
    B1 = DD;
end
SS=DD-B1;

