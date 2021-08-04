function [SS]=SSD_corr_rev4(CurrD,NextD,grid)

[l,c]=size(CurrD);
X=ones(l+(2*grid),c+(2*grid))*mean2(CurrD);
X(grid+1:l+grid,grid+1:c+grid)=CurrD;
SS=zeros((2*grid)+1,(2*grid)+1,l,c,'single');
for aa=1:(2*grid)+1
    for bb=1:(2*grid)+1
%         QR=X(aa:l+aa-1,bb:c+bb-1)-NextD;
%         QQ=QR.^2;
        SS(aa,bb,:,:)=movsum(movsum(((X(aa:l+aa-1,bb:c+bb-1)-NextD).^2),grid,1),grid,2);
%         SS(aa,bb,:,:)=S;
    end
end