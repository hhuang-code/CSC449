function v = hog(im,x,y,Wfull)

% function v = hog(im,x,y,Wfull)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
%  Compute the histogram of oriented gradidents on gray image (im) 
%   for a given location (x,y) and scale (Wfull)
%
%  v is the output column vector of the hog.  
%
%  Use Lowe IJCV 2004 Sections 5 and 6 to (1) adapt to local rotation
%    and (2) compute the histogram.  Use the parameters in the paper
%    Within the window a 4 by 4 array of histograms of oriented gradients 
%    with 8 discretized orientations per bin.  Do it separately per color channel
%    and then concatenate the resulting vectors.
%    Each v should be 3*128 dimensions = 3*4*4*8.
%    Finally, we will return a single concatenated vector: v = v(:).
%

v = zeros(3, 128); 

%%%%%%%%  fill in below

[x y] = deal(y, x);

[~, ~, chs] = size(im);  % number of channels

L = zeros(size(im)); % blurred image
for ch = 1 : chs
    L(:, :, ch) = conv2(im(:, :, ch), fspecial('gaussian', Wfull, 2.5), 'same');
end

% to avoid kepoints near image boarder, padding the image
L = padarray(L, [ceil(Wfull / 2) ceil(Wfull / 2)], 0, 'both');
[x y] = deal(x + ceil(Wfull / 2), y + ceil(Wfull / 2));

% calculate x derivative and y derivative
dx = zeros(size(L));
dy = zeros(size(L));
for ch = 1 : chs
    dx(:, :, ch) = conv2(L(:, :, ch), fspecial('sobel'), 'same');
    dy(:, :, ch) = conv2(L(:, :, ch), fspecial('sobel'), 'same');
end

% calculate magnitude and orientation
mag = zeros(size(L));
ori = zeros(size(L));
for ch = 1 : chs
    mag(:, :, ch) = sqrt(dx(:, :, ch) .* dx(:, :, ch) + dy(:, :, ch) .* dy(:, :, ch));
    ori(:, :, ch) = atan2(dy(:, :, ch), dx(:, :, ch)) + pi; % [0, 2 * pi]
end

% Gaussian window for weighted magnitude
win = fspecial('gaussian', Wfull, 1.5);

% top-left of the region
[x1 y1] = deal(x - ceil(Wfull / 2), y - ceil(Wfull / 2));
% bottom-right of the region
[x2 y2] = deal(x + ceil(Wfull / 2) - 1, y + ceil(Wfull / 2) - 1);

% orientation interval
interval = 2 * pi / 8;

% for each channel
for ch = 1 : chs
    tmp = mag(:, :, ch);
    rweight = tmp(x1 : x2, y1 : y2) .* win;
    tmp = ori(:, :, ch);
    rori = tmp(x1 : x2, y1 : y2);
    hist = zeros(4, 4, 8);  % eight bins for each blocks; 16 blocks
    % for each block
    for i = 1 : 4
        for j = 1 : 4
            % for each location in a block
            for r = (i - 1) * 4 + 1 : i * 4
                for c = (j - 1) * 4 + 1 : j * 4
                    if rori(r, c) == 2 * pi % the orientation of 2 * pi is the same as 0
                    rori(r, c) = 0;
                    end
                    hist(i, j, floor(rori(r, c) / interval) + 1) = ...
                        hist(i, j, floor(rori(r, c) / interval) + 1) + rweight(r, c);
                end
            end
            [~, idx] = max(hist(i, j, :));
            hist(i, j, :) = circshift(hist(i, j, :), -(idx - 1));
        end
    end
    v(ch, :) = hist(:);
end

%%%%%%%%  fill in above

v = v(:);