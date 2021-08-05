function [s,Pf]=powerSpectrum(Abkgs,pxlsize)
% computes 1D power spectrum of image Abkgs 
[N,M] = size(Abkgs);
Mn=mean2(Abkgs);
Img=Abkgs-Mn;
imgf = 10*log10(abs(fftshift(fft2(Img)).^2));
imgfp = ((imgf)/(N*M));

%%
% https://www.mathworks.com/matlabcentral/fileexchange/23636-radially-averaged-power-spectrum-of-2d-real-valued-matrix?focused=5188014&tab=function&s_tid=gn_loc_drop

dimDiff = abs(N-M);
dimMax = max(N,M);
if N > M                                                                    % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(N,dimDiff/2) imgfp NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        imgfp = [NaN(N,floor(dimDiff/2)) imgfp NaN(N,floor(dimDiff/2)+1)];
    end
elseif N < M                                                                % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        imgfp = [NaN(dimDiff/2,M); imgfp; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        imgfp = [NaN(floor(dimDiff/2),M); imgfp; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
end

halfDim = floor(dimMax/2) + 1;                                              % Only consider one half of spectrum (due to symmetry)
[X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
rho = round(rho);
i = cell(floor(dimMax/2) + 1, 1);
for r = 0:floor(dimMax/2)
    i{r + 1} = find(rho == r);
end
Pf = zeros(1, floor(dimMax/2)+1);
for r = 0:floor(dimMax/2)
    Pf(1, r + 1) = nanmean( imgfp( i{r+1} ) );
end
Fs_x=1/pxlsize; 
dx = 1/Fs_x;              % micrometers per pixel, pixel size in um
[M,N,~] = size(Img);      % pixels
x = dx*(0:N-1)';          % scale of image in um
dFx = Fs_x/N;             % cycles per micrometer, frequency increment
Fx = (-Fs_x/2:dFx:Fs_x/2-dFx)';     % cycles per micrometer, spatial frequency=2*pi/length
s=Fx(dimMax/2:dimMax);      % spatial frequancy distribution
end