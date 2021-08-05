%% MSD analysis from QPV velocity results
fdir='.\testdir\';   % data with phase images
fname=sprintf('Velws1_cell12.mat');  
load([fdir fname]);                 % load result from QPV velocity calc
sz=size(XS);
pxlsize=0.238;
for pp=1:sz(1)-1
    for qq=1:sz(2)-1
        x=(reshape(XS(pp,qq,:),1,sz(3))).*pxlsize;
        y=(reshape(YS(pp,qq,:),1,sz(3))).*pxlsize;
        t=Time;                                     % in minutes
        for kk=1:length(x)-1
            for jj = 1:length(x)-kk
               r2(jj) =((x(jj+kk)-x(jj)).^2 + (y(jj+kk)-y(jj)).^2);
               t2(jj) = t(jj+kk)-t(jj);
            end
            R2(kk)=mean(r2);    
            T2(kk)=mean(t2);
            clear t2 r2;
        end
        Diff{pp,qq}=R2;
        Timed{pp,qq}=T2;
        fitData=polyfit(log(T2),log(R2),1);
        if ~isnan(fitData(1)) && ~isnan(fitData(2))
            DiffC(pp,qq)=fitData(1);                % diffusion coefficient
            Expo(pp,qq)=(exp(fitData(2)))/4;        % anomalous constant
        else
            DiffC(pp,qq)=0;
            Expo(pp,qq)=0;
        end
    end
end

fname2=sprintf('MSDws12_1.mat');
save([fdir1 fname2],'Abkg_stored2','Diff','DiffC','Time','Timed','Expo');

%% MSD analysis results for diffusion coefficient and cell viscosity measurement
% load([fdir fname2],'DiffC','Abkg_stored2','Expo','Diff','Timed')
DD=Abkg_stored2(:,:,1);
sz1=size(Abkg_stored2);
SSf = imfilter(abs(DD), fspecial('gaussian', [20 20], 2));
[junk threshold] = edge(SSf, 'sobel');
fudgeFactor = 0.8;          % change to fit cell type
BWs = edge(SSf,'sobel', threshold * fudgeFactor);
se90 = strel('line', 15, 90);
se0 = strel('line', 15, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');
seD = strel('diamond',5);
BWfinal = imerode(BWdfill,seD);
BWfinal=imclearborder(BWfinal,4);
BWfinalCyto = bwareaopen(BWfinal, 2000);  % cell mask

aa=0; 
for pp=1:sz1(1)-1
    for qq=1:sz1(2)-1
        if BWfinalCyto(pp,qq)==1
            aa=aa+1;
            MSDCell(aa,:)=Diff{pp,qq};  % accumulate interrogation window MSD vs time lag 
        end
    end
end
MSDCellAvg=mean(MSDCell,1);             % average MSD from all windows
for dd=1:length(MSDCellAvg)
    MassTAll(dd,1)=mean(nonzeros(MSDCellAvg(:,dd)));
    SemErrAll(dd,1)=std(nonzeros(MSDCellAvg(:,dd)))/sqrt(nnz(MSDCellAvg(:,dd)));
end
TimeSec=Timed{pp,qq};                   % time vector same for all windows
fitData=polyfit(log(TimeSec(1:20)),log(MassTAll(1:20))',1);
AnoCell=fitData(1);                     % average diffusion coefficient from all windows
DiffCell=(exp(fitData(2)))/4;           % anomalous constant from all windows

[s,Pf]=powerSpectrum(Abkg_stored2(:,:,1),pxlsize); % power spectrum density from QPI image
LngS=sum(s>0.14); EffFq=sum((s(257-LngS+1:257)'.*Pf(257-LngS+1:257)))/sum(Pf(257-LngS+1:257));
EffDp=1/(EffFq);                        % effective particle size from power spectrum
k=1.38*(10^-23);                        % Boltzmann constant
T=310;                                  % temp in Kelvin
mu=(k*T*(10^18))./(3*3.14*EffDp*(DiffCell/60)); % cell viscosity
