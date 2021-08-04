function [BWfinal,SS] = imagebackground_poly4(I)
%function [B] = imagebackground_poly4(I)
%function to find the background of an image, I, using a 4th order
%polynomial fit to the background pixels
%input: I, the grayscale image to find the background of
%output: B, the background of I
%method: find 'objects' in I, mask them from the image, paint the remaining
%area using the inpaint_nans function

%step 2, detect entire cell
% Io=imtophat(single(I),strel('sphere',10));
% I = imfilter(abs(I), fspecial('gaussian', [10 10], 2));
% If1 = imfilter(I, fspecial('sobel'));
% If2 = imfilter(I, fspecial('sobel')');
% BWs = adaptivethreshold(abs(If1),[300 300], 0) | adaptivethreshold(abs(If2),[300 300], 0);

[junk threshold] = edge(I, 'sobel');
fudgeFactor = 0.8; %was 0.4 for RBCs
BWs = edge(I,'sobel', threshold * fudgeFactor);
% BWs = adaptivethreshold(abs(I),[100 100], 0);
% figure(2)
% imshow(BWs)
% title('binary gradient mask');

%step 3, dilate the image
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);

BWsdil = imdilate(BWs, [se90 se0]);
% figure(3)
% imshow(BWsdil)
% title('dilated gradient mask')

%step 4, fill gaps
BWdfill = imfill(BWsdil, 'holes');
% figure(4)
% imshow(BWdfill);
% title('binary image with filled holes');



%step 6, smooth image
seD = strel('diamond',5);
%BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 2000); %remove regions smaller than 10 pixels
se90 = strel('line', 18, 90);
se0 = strel('line', 18, 0);
BWfinal=imdilate(BWfinal,[se0 se90]);

IList = I(~BWfinal);

sz = size(I);
[X,Y] = meshgrid(1:sz(2), 1:sz(1));
XList = X(~BWfinal);
YList = Y(~BWfinal);

if sum(~isnan(BWfinal)) ~=0
    CFit = polyfitn([XList,YList], IList, 4);
    B =(reshape(polyvaln(CFit, [X(:), Y(:)]), sz(1), sz(2)));
else
    B = I;
end
SS=B-I;
% M = ~BWfinal; %return the mask used for processing
