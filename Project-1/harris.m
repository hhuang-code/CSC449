function C = harris(dx,dy,Whalfsize)

% function C = harris(dx,dy,Whalfsize)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
%   dx is the horizontal gradient image
%   dy is the vertical gradient image
%   Whalfsize is the half size of the window.  Wsize = 2*Whalfsize+1
%
%  output
%   C is an image (same size as dx and dy) where every pixel contains the
%   corner strength.  



% kappa = ones(2*Whalfsize+1);

%%%%%%%%% fill in below

win = fspecial('gaussian', 2 * Whalfsize + 1, 0.7);

m11 = conv2(dx .* dx, win, 'same');
m12 = conv2(dx .* dy, win, 'same');
m21 = conv2(dx .* dy, win, 'same');
m22 = conv2(dy .* dy, win, 'same');

[row, col] = size(m11);

C = ones(size(m11));
% Shi-Tomasi method
for i = 1 : row
    for j = 1 : col
        M = [m11(i, j) m12(i, j); m21(i, j) m22(i, j)];
        C(i, j) = min(eig(M));
    end
end

%%%%%%%% done
