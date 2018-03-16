function [X,Sout] = onedetect(index,Whalfsize,method)

% function onedetect(index,Whalfsize)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Wrapper for function to detect interest points in the image
%
%
%  Will run harris or dog detector depending on arguments.
%    run as X = onedetect(index,Whalf) will just run harris.
%    run as [X,S] = onedetect(index,Whalf,'dog') will run dog
%      S is the output scale of a point
%
%  index is the index for the image
%  Whalfsize is the half size of the window.  Wsize = 2*Whalfsize+1
%
%  output
%  X is a 2xn matrix of corners locations where n is induced by the data
%         top row is the horizontal,x coordinate (column, not row)
%  Sout (optional) is the scale of the detected points if the method can do it

global cards;

if nargin == 1
    Whalfsize=7;
    method = 'harris';
end

im = cards{index};
gm = rgb2gray(im);

if strcmp(method,'harris')
 
    dy = conv2(gm,fspecial('sobel'),'same');
    dx = conv2(gm,fspecial('sobel')','same');
    
    C = harris(dx,dy,Whalfsize);
    
    [Xr,Xc] = nonmaxsuppts(C,3,0.2);
    
    %prune border
    [r,c] = size(gm);
    by = Xr > 16 & Xr < r-16;
    bx = Xc > 16 & Xc < c-16;
    Xr = Xr(by & bx);
    Xc = Xc(by & bx);
    
    X = [Xc';Xr'];

elseif strcmp(method,'dog')
    
    [X,S] = dog(gm, 16);
    
    if (nargout == 2)
        Sout = S;
    end
end

