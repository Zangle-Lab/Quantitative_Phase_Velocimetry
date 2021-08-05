%%  Measurement of intracellular mass transport velocity using quantitative phase velocimetry
% main script for performing QPV on QPI images
% SP 13th June 2019
% first section: define image processing parameters (file locations, 
% processing options, etc.)
% second section: load phase images, the time captured
% third section: initialie grid parameters overlayed on image
% fourth section: background correct phase images to flaten
% fifth section: perform QPV using SSD on images, store displacements and
% deform initialized grid with the displacements
% sixth section: display the velocity on images over time

clc; 
close all;
clear all;

%% section 1: set the parameters for computation
fdir='.\testdir\';   % data with phase images
xcg=4;               % window size for spatial averaging of SSD
gs=15;               % SSD interrogation window size, chose an odd number if possible
sgap=4;              % gap bwn frames undergoing SSD, must be less than numf
numf=20;             % total number of frames considered for calculation
medfiltsz=4*xcg;     % size of selective median filter
pxl_conv=0.238;      % width of image pixel for 120X
Ref_inc=0.18;        % refractive index increment in um3/pg
oplf=0.623;          % optical path length factor to convert QPI phase images to dry mass distribution
Pixel_area=pxl_conv*pxl_conv; % pixel area of image
fstart=1;            % First frame number
mm=12;               % imaging position number

%% section 2: stored the phase images and their time
% time stamp is taken from tif image as phase image was produced later from tif
% images using Phasics SDK library
fname=sprintf('Cell_%d_pos_%d.mat',mm,fstart);   % first phase image for initialization            
fname1=sprintf('Cell_%d_pos_%d.tif',mm,fstart);  % tif image for time stamp
D=dir([fdir fname1]);           
Timenum(1)=D.datenum;                            % time stamp for first image frame
lnumf=ceil(numf./sgap);                          % total number of frames processed
load([fdir fname],'Phase');
for uu=1:ceil(numf./sgap)
    fname=sprintf('Cell_%d_pos_%d.mat',mm,(((uu-1)*sgap)+fstart));  
    load([fdir fname],'Phase');
    fname1=sprintf('Cell_%d_pos_%d.tif',mm,(((uu-1)*sgap)+fstart));
    D=dir([fdir fname1]);
    Timenum(uu+1)=D.datenum;
    Time(uu)=(Timenum(uu+1)-Timenum(1))*24*60;   %time in min
    D_stored2(:,:,uu)=Phase;                     % 3D array with phase images over time
end
clear Abkg Phase L R B M;

%% section 3: initialise the computation by defining the zero vectors to which the results are stored
sz = size(D_stored2);
kkf = sz(3)-1;                            % index of last frame to be used as CurrD image in velocity computation
[X0,Y0]=meshgrid(1:sz(1)+1,1:sz(1)+1);
sX = single(X0);
sY = single(Y0);
XS = single(X0);                          % initialize x grid uniform grid at time=0 min
YS = single(Y0);                          % initialize y grid uniform grid at time=0 min
dX=zeros((sz(1)-gs),(sz(2)-gs),'single'); % initialize x dispalcement matrix
dY=dX;                                    % initialize y displacemenet matrix
clear X0 Y0 sX0 sY0; 
    
%% section 4: Background correction of QPI phase images
A_stored2=zeros(sz,'single');             
Abkg_stored2=A_stored2;                   % initialize 3D array for storing backround corrected images
for ss=1:lnumf
    [B,A_stored2(:,:,ss)]=imagebackground_poly4(D_stored2(:,:,ss));  % background correcting phase image
    Abkg_stored2(:,:,ss)= A_stored2(:,:,ss)*(oplf*Pixel_area/Ref_inc); % convert phase image to dry mass map    
    mass(ss)=sum(sum(Abkg_stored2(:,:,ss).*B));                       % total dry mass of cells in frame
end

%% section 5: perform QPV on image frames and use displacement result to deform grid from t=0
for kk=1:kkf 
      CurrD = gpuArray(Abkg_stored2(:,:,kk));   % SSD is perform between CurrD and NextD
      NextD = gpuArray(Abkg_stored2(:,:,kk+1)); % remove gpuArray conversion if not supported
      XCs=SSD_corr_rev9(NextD,CurrD,gs);        % SSD computation, use SSD_corr_rev4 for non-GPU SSD computation
      clear CurrD NextD;
      XCs=gather(XCs);
      XCsp=zeros((2*gs)+1,(2*gs)+1,sz(1),sz(2),'single');
      for uu=1:sz(1)-xcg+1    %spatial averaging of SSD
          for vv=1:sz(2)-xcg+1
              XCsp(:,:,uu,vv)=mean(mean(XCs(:,:,uu:uu+xcg-1,vv:vv+xcg-1),3),4);  % averaging xcg by xcg neighborhood SSD planes to reduce noise
          end
      end
      clear XCs;
      dX1=zeros(sz(1)+1);
      dY1=dX1;
      for ii =1:sz(1)-3
           for jj =1:sz(2)-3
              [xm,ym]=findvalley_v3_rev4(XCsp(:,:,ii,jj));  % find lowest pixel & lowest point inside the pixel
              dX1(ii+2,jj+2) = xm;                          % x displacement in pixels
              dY1(ii+2,jj+2) = ym;                          % y displacement in pixels
           end
      end
    clear XCsp; 
    % use selective median filter to remove spurious velocities and not
    % alter the original values
    for ee=medfiltsz+1:sz(1)-medfiltsz+1
        for ff=medfiltsz+1:sz(2)-medfiltsz+1
            Ux=dX1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
            ux=median(Ux(:),'omitnan');
            Uy=dY1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
            uy=median(Uy(:),'omitnan'); 
            if abs(dX1(ee,ff))>(1.5*ux) || abs(dX1(ee,ff))<(0.7*ux)
                dX(ee,ff,kk-tcg+1)=ux;
            else
                dX(ee,ff,kk-tcg+1)=dX1(ee,ff);
            end
            if abs(dY1(ee,ff))>(1.5*uy) || abs(dY1(ee,ff))<(0.7*uy)
                dY(ee,ff,kk-tcg+1)=uy;
            else
                dY(ee,ff,kk-tcg+1)=dY1(ee,ff);
            end
        end
    end
    for pp=medfiltsz+1:sz(1)-medfiltsz+1
        for qq=medfiltsz+1:sz(2)-medfiltsz+1
            sx=sX(pp,qq);
            rx=floor(sx);
            sy=sY(pp,qq);
            ry=floor(sy);
            if rx>0 && rx<sz(1)-medfiltsz
                if ry>0 && ry<sz(1)-medfiltsz
                    sX(pp,qq) = sX(pp,qq) + (dX(ry,rx,kk));  % add the disaplacement to existing grid to deform based on cell movement
                    sY(pp,qq) = sY(pp,qq) + (dY(ry,rx,kk));
                end
            end
        end
    end
    XS(:,:,kk+1) = sX;          % Store the deformed grid over time
    YS(:,:,kk+1) = sY;
    fprintf('Completed loop %d\n',kk);
    clear CurrD NextD dX1 dY1;
end
fnames=sprintf('Velws1_cell12.mat');
save([fdir fnames]);

%% section 6: display the velocity vectors calculated
figure(1);
se90 = strel('line', 35, 90);
se0 = strel('line', 35, 0); 
for zz=1:kkf
    [BWfinal,A_stored2(:,:,zz)]=imagebackground_poly4(D_stored2(:,:,zz));
    BWmap=imerode(BWfinal(8:504,8:504,zz),[se0 se90]);
    DX=imresize((BWmap.*dX(:,:,zz)),0.124);                     % mask out the vectors outside the cell
    DY=imresize((BWmap.*dY(:,:,zz)),0.124);
    x=[(-sz(2)/4)+4 , (sz(2)/4)-3];                                                      % x axis ticks with origin as center
    y=[(-sz(1)/4)+4 , (sz(1)/4)-3]; 
    [Vxo,Vyo]=meshgrid((-sz(2)/4)+4:xcg:(sz(2)/4)-4-1,(-sz(1)/4)+4:xcg:(sz(1)/4)-4-1);   % center axis at origin to quiver
    imagesc(x,y,Abkg_stored2(:,:,zz)); 
    hold on;  caxis([-0.01 0.05]);
    hold on 
    DX(57,10)=1/pxl_conv;                     % 1 um scale
    h1=quiver(Vxo,Vyo,DX,DY,0,'k'); axis off; % turn scale off to get absolute value of velocity vectors
    scale=3;
    hU1 = get(h1,'UData');
    hV1 = get(h1,'VData');
    set(h1,'UData',scale*hU1,'VData',scale*hV1)
    h1.LineWidth=1;
    textside=sprintf('1 µm/min',zz);
    text('units','pixels','position',[50 110],'fontsize',25,'color','w','string',textside)
    hold off
    set(gcf,'color','w'); pause(1);
end