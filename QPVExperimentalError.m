%% Experimental error chart modelling : SP 
% set the parameters for computation
    xcg=4;              % window size for spatial averaging of SSD
    gs=3:2:57;          % SSD window size, chose an odd number
    w=2*gs;             % SSD search window size in pixel
    gap=3:2:401;        % difference bwn frames undergoing SSD, must be less than (numf-tcg)
WindowSize=w.*0.238; 	% converted to um
VelocityY=(gap-1).*0.05; % in um, cell moved vertically 0.05 um every frame
VelocityX=(gap-1).*0;    % in um
medfiltsz=4;             % selective median filter size

%% initialise the computation by defining the zero vectors to which the results are stored
fdir='.\testdir\';
load('.\testdir\Cell_12_pos_1.mat');        % first cell image
[B,massA]=imagebackground_poly4(Phase);     % image background correction
CurrD=massA;
sz = size(massA);
MassA=abs(imresize(massA((2*xcg)+1:sz(1)-(2*xcg),(2*xcg)+1:sz(1)-(2*xcg)),(1/xcg)));
MA=abs(sum(sum(MassA)));
c = kmeans(massA(:), 2, 'MaxIter', 10000);
CC=reshape(c,sz(1),sz(2));
BW=CC~=mode(CC);
BWmap=abs(imresize(BW((2*xcg)+1:sz(1)-(2*xcg),(2*xcg)+1:sz(1)-(2*xcg)),(1/xcg))); % cell mask

%% compute the displacement using QPV and then dispalcemeent error
for tt=1:length(gs)
    for kk=1:length(gap)
        fname=sprintf('Cell_12_pos_%d.mat',gap(1,kk));
        load([fdir fname],'Phase');
        [B2,massB]=imagebackground_poly4(Phase);
        NextD=massB;
        XCs=SSD_corr_rev4(NextD,CurrD,gs(1,tt));
        XCsp=zeros((2*gs(1,tt))+1,(2*gs(1,tt))+1,sz(1)/xcg,sz(2)/xcg);
        for uu=1:xcg:sz(1)-xcg+1                 %spatial averaging of SSD
            for vv=1:xcg:sz(2)-xcg+1
                XCsp(:,:,((uu-1)/xcg)+1,((vv-1)/xcg)+1)=mean(mean(XCs(:,:,uu:uu+xcg-1,vv:vv+xcg-1),3),4);
            end
        end
        clear NextD XCs Phase Abkg fname uu vv;
        for ii =1:floor(sz(1)/xcg)              % finding lowest pixel
            for jj =1:floor(sz(2)/xcg)
                [dX1(ii,jj),dY1(ii,jj)] = findvalley_v3_rev4(XCsp(:,:,ii,jj));
            end
        end
        clear XCsp;
        % selective median filtering to remove spurious velocities
        for ee=medfiltsz+1:floor(sz(1)/xcg)-medfiltsz
            for ff=medfiltsz+1:floor(sz(2)/xcg)-medfiltsz
                Ux=dX1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
                ux=median(Ux(:),'omitnan');
                Uy=dY1(ee-medfiltsz:ee+medfiltsz,ff-medfiltsz:ff+medfiltsz);
                uy=median(Uy(:),'omitnan'); 
                if abs(dX1(ee,ff))>(1.5*ux) || abs(dX1(ee,ff))<(0.7*ux)
                    dX(ee,ff)=ux;
                else
                    dX(ee,ff)=dX1(ee,ff);
                end
                if abs(dY1(ee,ff))>(1.5*uy) || abs(dY1(ee,ff))<(0.7*uy)
                    dY(ee,ff)=uy;
                else
                    dY(ee,ff)=dY1(ee,ff);
                end
            end
        end
        
    %% Displacemnt magntiude and direction error compared to actual displacement
        if sum(sum(BWmap.*dX))~=0
            Vxavg(kk,tt)=mean2(nonzeros(dX.*BWmap)); % error in x  
        else
            Vxavg(kk,tt)=0;
        end
        if sum(sum(BWmap.*dY))~=0
            Vyavg(kk,tt)=mean2(nonzeros(BWmap.*dY));  % error in y
        else
            Vyavg(kk,tt)=0;
        end
        T(kk,tt)=(atan2d(Vyavg(kk,tt),Vxavg(kk,tt)));                   % dispalcement direction
        Vone(kk,tt)=(sqrt((Vxavg(kk,tt)^2)+(Vyavg(kk,tt)^2)))*0.238;    % displacemenet magntide from x and y components   
        Mag_accuracy(kk,tt)=((VelocityY(1,kk)-Vone(kk,tt))/VelocityY(1,kk))*100;    % magnitude error
        Dir_accuracy(kk,tt)=(90-T(kk,tt));                              % direction error

    %% save error results
        fprintf('Completed loop %d, %d \n',kk,tt);       
        fname1=sprintf('ExpVelocityError.mat');
        save([fdir fname1]);
    end
end

%% Display magntide and direction error data
Window=(3:2:57).*0.238*2;
DisplacementRange=(gap-1)*0.05;
figure(1)
imagesc(Window, DisplacementRange,abs(Mag_accuracy)); % magnitude error chart
colormap(flipud(hot))
caxis([0 50]);
set (gca,'Ydir','normal')
xlabel('Window size (\mum)');
ylabel('Displacement (\mum)');
colorbar;
figure(2)
imagesc(Window, DisplacementRange,abs(Dir_accuracy))  %direction error chart
colormap(flipud(hot))
caxis([0 40]);
set (gca,'Ydir','normal')
xlabel('Window size (\mum)');
ylabel('Displacement (\mum)');
colorbar;