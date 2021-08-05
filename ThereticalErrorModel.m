%% Model for predicting displacement magntiude and direction error from QPV
% load cell QPI image
load('./testdir/Cell_12_pos_1.mat');
[B,~]=imagebackground_poly4(Phase); % background correction
Abkg1=B-Phase;
[Abkg2,BWmap]=imagebackground_poly4_kmeans(Abkg1); % further background correction if required
Abkg=Abkg2;
pxlsize=0.238; % pixel width of image in um

%% mask background to get power spectrum of just cell
BWsdil=Abkg>0.005;
BWdfill = imfill(BWsdil, 'holes');
seD = strel('diamond',2);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000);
se90 = strel('line', 5, 90);
se0 = strel('line', 5, 0);
Maskw = imdilate(BWfinal, [se90 se0]);
Abkgs=Maskw.*Abkg;

%% power spectrum image of cell
[s,Pf]=powerSpectrum(Abkgs,pxlsize);
S=1./s;                   %length (um), in 2pi, 2times diameter of particles

%% use power spectral density to find fraction of particles of different diameter
AreaPSD5=sum(Pf(5:end));
FractionPSD5=Pf(5:end)./AreaPSD5;
Pf5=Pf;
FractionPSD6=flip(FractionPSD5);
Fs_x=1/pxlsize;
dx = 1/Fs_x;              % micrometers per pixel, pixel size in um

%% Point spread function of the optics used
M=120;                                  % objective magnification
NA=1.3;                                 % numerical aperture of objective
lamda=0.623;                            % illumination wavelength in um
Window=(3:2:57).*dx.*2;
dsObj=((1.22*(M+1)*lamda)/(NA))/2;      % PSF diameter by 2 as it is phase image in um
DiffLim=lamda/NA;                       % diffraction limit on object plane

%% comparison to experimental results
load('./testdir/ExpVelocityError.mat');         % load experiemntal error data
gap=2:2:200;
gs=3:2:143;
gsum=gs*0.238;
WindowRange=Window;
DisplacementRange=(gap-1)*0.05;

%% Construction of magnitude error chart based on above values
DisRange=DisplacementRange;
W0=zeros(length(DisRange),length(WindowRange),length(S)-4);    %initialize matrix for storing error values
W1=zeros(length(DisRange),length(WindowRange));
for pp=1:length(S)-4
    ParticleDia(pp)=S(length(S)-pp+1);                         % particle size obtained from power spectrum in 2pi, dia is pi
    for dd=1:length(DisRange)
        Displacement=DisRange(dd);
        for ww=1:length(WindowRange)
            WinSize=WindowRange(ww);
            if Displacement>=(WinSize/2)-(ParticleDia(pp)/2)   % if displacement cannot be seen by interrogation window used
                W1(dd,ww)=(100-abs(FractionPSD6(pp)*100));
            end
            if Displacement<=(WinSize/2)-(ParticleDia(pp)/2)   % if dispalcement is visible by window
                W1(dd,ww)=abs(FractionPSD6(pp)*100);
            end
            if Displacement<=DiffLim                           % if displacement below diffraction limit of optics
                W1(dd,ww)=100;
            end
        end
    end
    W0(:,:,pp)=W1;
end
WWW=zeros(length(DisRange),length(WindowRange));               % average over all particle sizes keeping particle sizes visible by window size
for ww=1:length(WindowRange)
    QW=ParticleDia<WindowRange(ww); PDRange=sum(QW);
    WWW(:,ww)=mean(W0(:,ww,1:PDRange),3);
end

%% display comparison between experiemntal and theoretical direction error
% theoretical magnitude error display
figure(2);
imagesc(WindowRange,DisRange,abs(WWW));
colormap(flipud(hot))
caxis([0 100]);
set (gca,'Ydir','normal')
xlabel('Window size (\mum)');
ylabel('Displacement (\mum)');
colorbar;
title('Magnitude Theoretical error');
% experimental magnitude error display
figure(3);
imagesc(WindowRange,DisplacementRange,abs(Mag_accuracy(1:100,:)));
colormap(flipud(hot))
caxis([0 100]);
set (gca,'Ydir','normal')
xlabel('Window size (\mum)');
ylabel('Displacement (\mum)');
colorbar;
title('Magnitude Experimental error');
set(gcf,'color','w');

%% Direction error chart
Wd=zeros(length(DisRange),length(Window),length(S)-2);    %initialize matrix for storing error values
W2=zeros(length(DisRange),length(Window));
for aa=1:length(S)-2
    ParticleDia(aa)=S(length(S)-aa+1);
    for bb=1:length(DisRange)
        Displacement=DisRange(bb);
        for cc=1:length(Window)
            WinSize=Window(cc);
            if (WinSize/2)<Displacement-ParticleDia(aa) || (WinSize/2)<ParticleDia(aa)/sqrt(2) % displacement of particle not visible by used window size
                W2(bb,cc)=100;
            end
            if (WinSize/2)>Displacement-ParticleDia(aa) && (WinSize/2)>ParticleDia(aa)/sqrt(2) % if displacement visible by window used
                W2(bb,cc)=0;
            end
            if Displacement<=DiffLim % if displacement below diffraction limit
                W2(bb,cc)=100;
            end
        end
    end
    Wd(:,:,aa)=W2;
end
WWd=zeros(length(DisRange),length(WindowRange));     % average over all particle sizes keeping particle sizes visible by window size
for ww=1:length(WindowRange)
    QW=ParticleDia<WindowRange(ww); PDRange=sum(QW);
    WWd(:,ww)=mean(Wd(:,ww,1:PDRange),3);
end

%% display experimental and theoretical direction error
figure(5);
imagesc(Window,DisRange,WWd);
colormap(flipud(hot))
caxis([0 50]);
set (gca,'Ydir','normal')
xlabel('Window size (\mum)');
ylabel('Displacement (\mum)');
colorbar;
title('Direction Theoretical error');
set(gcf,'color','w');

figure(6);
imagesc(WindowRange,DisplacementRange,abs(Dir_accuracy(1:100,:)));
colormap(flipud(hot))
caxis([0 50]);
set (gca,'Ydir','normal')
xlabel('Window size (\mum)');
ylabel('Displacement (\mum)');
colorbar;
title('Direction Experimental error');
set(gcf,'color','w');
