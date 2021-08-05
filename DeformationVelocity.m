%% deformation velocity from QPV velocities
fdir='.\testdir\'; 
fname=sprintf('Velws1_Cell12');
load([fdir fname]);                         % load results from QPV main code
pxlsize=0.238;                              % pixel width for our 120X mag
DD=Abkg_stored2(:,:,1);
c = kmeans(DD(:), 3, 'MaxIter', 10000);     % cell segmentation
CC=reshape(c,512,512);
BWs1=CC~=mode(CC(:));
BWdfill = imfill(BWs1,8, 'holes');
seD = strel('diamond',1);
%BWfinal = imerode(BWnobord,seD);
BWfinal = imdilate(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000);
se90 = strel('line', 12, 90);
se0 = strel('line', 12, 0);
BWfinal=imerode(BWfinal,[se0 se90]);        % constrict mask
 BWfinal = bwareaopen(BWfinal, 5000);
 se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);       % final cell mask

Vel=(sqrt((dX(:,:,1)).^2+(dY(:,:,1)).^2))*pxlsize;      % intracellular cell  magnitude map 
VelCell=mean2(nonzeros(Vel.*BWfinal(8:504,8:504)));     % whole cell velocity 
DefVel=(abs((Vel-VelCell)));                            % Cell deformation velocity map