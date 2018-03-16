function [X, S] = dog(gm, Wfull)

% function [X, S] = dog(gm)
%
%   gm is an input grayscale image
%   Wfull is scale, it should be 16 * (2 ^ k)
%
%   output
%       X is a 2xn matrix of corners locations where n is induced by the data
%           top row is the horizontal,x coordinate (column, not row)
%       Sout (optional) is the scale of the detected points if the method can do it
%


assert(floor(log2(Wfull / 16)) == log2(Wfull / 16));

target_octave = floor(log2(Wfull / 16)) + 1;

% define parameters
num_octaves = 3; % number of octaves
s = 2;  % called intervals, produce (s + 3) blurred images per octave_config_info
k = 2 ^ (1.0 / s);  % multiplicative factor
init_sigma = 1.6;  % initial sigma for each octave
init_win = fspecial('gaussian', ceil(6 * init_sigma), init_sigma);  % initial guassian
threshold = 0.03;   % lower constrast keypoint will be discarded
r = 10; % curvature threshold

% no more than 3 octaves
assert(target_octave <= num_octaves);

% allocate memory for blurred images and dog images
gaussians = cell(num_octaves, s + 3); % each row for an octave
dogs = cell(num_octaves, s + 2);  % each row for an octave

% bottom image in the first octave
gaussians(1, 1) = {(conv2(gm, init_win, 'same'))};

for i = 1 : num_octaves
    sigma = init_sigma; % reset gaussian sigma
    for j = 2 : (s + 3)
        sigma = sigma * (k ^ (j - 1));
        win = fspecial('gaussian', ceil(6 * init_sigma), sigma);  % win size unchanged
        prev_image = cell2mat(gaussians(i, (j - 1)));  % previous blurred image
        curr_image = conv2(prev_image, win, 'same');  % current blurred image
        dog_image = curr_image - prev_image;
        % store current blurred image and dog image
        gaussians(i, j) = {curr_image};
        dogs(i, j - 1) = {dog_image};    
    end
    
    % construct the bottom image in the next octave
    if (i < num_octaves)
        gaussians(i + 1, 1) = {impyramid(cell2mat(gaussians(i, s + 1)), 'reduce')};
    end
end

X = []; % keypoints
for i = 1 : num_octaves
    octave_kp = [];  % keypoints detected in all dog images at i-th octave
    for j = 2 : (s + 3 - 1 - 1) % omit the top and bottom dog images
        curr_dog = cell2mat(dogs(i, j));
        below_dog = cell2mat(dogs(i, j - 1));
        up_dog = cell2mat(dogs(i, j + 1));
    
        [nrow ncol] = size(curr_dog); % dog image size
    
        % detect extrama excluding the dog image boarder
        area = curr_dog(2 : nrow - 1, 2 : ncol - 1);
    
        % -----------------detect local maxima ------------------
        top_left = area > curr_dog(1 : nrow - 2, 1 : ncol - 2);
        top = area > curr_dog(1 : nrow - 2, 2 : ncol - 1);
        top_right = area > curr_dog(1 : nrow - 2, 3 : ncol);
        left = area > curr_dog(2 : nrow - 1, 1 : ncol - 2);
        right = area > curr_dog(2 : nrow - 1, 3 : ncol);
        bottom_left = area > curr_dog(3 : nrow, 1 : ncol - 2);
        bottom = area > curr_dog(3 : nrow, 2 : ncol - 1);
        bottom_right = area > curr_dog(3 : nrow, 3 : ncol);
        % maxima compared with current image
        local_max = top_left & top & top_right & left & right & bottom_left & bottom & bottom_right;
    
        top_left = area > up_dog(1 : nrow - 2, 1 : ncol - 2);
        top = area > up_dog(1 : nrow - 2, 2 : ncol - 1);
        top_right = area > up_dog(1 : nrow - 2, 3 : ncol);
        left = area > up_dog(2 : nrow - 1, 1 : ncol - 2);
        middle = area > up_dog(2 : nrow - 1, 2 : ncol - 1);
        right = area > up_dog(2 : nrow - 1, 3 : ncol);
        bottom_left = area > up_dog(3 : nrow, 1 : ncol - 2);
        bottom = area > up_dog(3 : nrow, 2 : ncol - 1);
        bottom_right = area > up_dog(3 : nrow, 3 : ncol);
        % maxima compared with up image
        local_max = local_max & top_left & top & top_right & left & middle & right ... 
                & bottom_left & bottom & bottom_right;
                
        top_left = area > below_dog(1 : nrow - 2, 1 : ncol - 2);
        top = area > below_dog(1 : nrow - 2, 2 : ncol - 1);
        top_right = area > below_dog(1 : nrow - 2, 3 : ncol);
        left = area > below_dog(2 : nrow - 1, 1 : ncol - 2);
        middle = area > below_dog(2 : nrow - 1, 2 : ncol - 1);
        right = area > below_dog(2 : nrow - 1, 3 : ncol);
        bottom_left = area > below_dog(3 : nrow, 1 : ncol - 2);
        bottom = area > below_dog(3 : nrow, 2 : ncol - 1);
        bottom_right = area > below_dog(3 : nrow, 3 : ncol);
        % maxima compared with below image
        local_max = local_max & top_left & top & top_right & left & middle & right ...
                & bottom_left & bottom & bottom_right;
    
        % -----------------detect local minima ------------------
        top_left = area < curr_dog(1 : nrow - 2, 1 : ncol - 2);
        top = area < curr_dog(1 : nrow - 2, 2 : ncol - 1);
        top_right = area > curr_dog(1 : nrow - 2, 3 : ncol);
        left = area > curr_dog(2 : nrow - 1, 1 : ncol - 2);
        right = area > curr_dog(2 : nrow - 1, 3 : ncol);
        bottom_left = area > curr_dog(3 : nrow, 1 : ncol - 2);
        bottom = area > curr_dog(3 : nrow, 2 : ncol - 1);
        bottom_right = area > curr_dog(3 : nrow, 3 : ncol);
        % minima compared with current image
        local_min = top_left & top & top_right & left & right & bottom_left & bottom & bottom_right;
    
        top_left = area < up_dog(1 : nrow - 2, 1 : ncol - 2);
        top = area < up_dog(1 : nrow - 2, 2 : ncol - 1);
        top_right = area < up_dog(1 : nrow - 2, 3 : ncol);
        left = area < up_dog(2 : nrow - 1, 1 : ncol - 2);
        middle = area < up_dog(2 : nrow - 1, 2 : ncol - 1);
        right = area < up_dog(2 : nrow - 1, 3 : ncol);
        bottom_left = area < up_dog(3 : nrow, 1 : ncol - 2);
        bottom = area < up_dog(3 : nrow, 2 : ncol - 1);
        bottom_right = area < up_dog(3 : nrow, 3 : ncol);
        % minima compared with up image
        local_min = local_min & top_left & top & top_right & left & middle & right ...
                & bottom_left & bottom & bottom_right;
                
        top_left = area < below_dog(1 : nrow - 2, 1 : ncol - 2);
        top = area < below_dog(1 : nrow - 2, 2 : ncol - 1);
        top_right = area < below_dog(1 : nrow - 2, 3 : ncol);
        left = area < below_dog(2 : nrow - 1, 1 : ncol - 2);
        middle = area < below_dog(2 : nrow - 1, 2 : ncol - 1);
        right = area < below_dog(2 : nrow - 1, 3 : ncol);
        bottom_left = area < below_dog(3 : nrow, 1 : ncol - 2);
        bottom = area < below_dog(3 : nrow, 2 : ncol - 1);
        bottom_right = area < below_dog(3 : nrow, 3 : ncol);
        % minima compared with below image
        local_min = local_min & top_left & top & top_right & left & middle & right ...
                & bottom_left & bottom & bottom_right;
                
        % local extrama is either local maxima or local minima
        l_e = local_max | local_min;
    
        [x, y] = find(l_e); % coordinates of extramas in the 'area' of current dog image
        num_e = size(x);    % number of extramas in the defined 'area' of current dog image
    
        if (num_e >= 1) % current dog image has local extramas
            % Hesian matrix of current dog image
            [dx, dy] = gradient(double(curr_dog));
            [dxx, dxy] = gradient(dx);
            [dxy, dyy] = gradient(dy);
        end
    
        for p = 1 : num_e
            ex = x(p) + 1;  % x corrdinates of extramas in current dog image
            ey = y(p) + 1;  % y corrdinates of extramas in current dog image
            % discard lower constrast keypoint
            if (abs(curr_dog(ex, ey)) < threshold)
                l_e(ex - 1, ey - 1) = 0;    % pay attention to '-1'
            % eliminate edge responses
            else
                trace = dxx(ex, ey) + dyy(ex, ey);
                determinant = dxx(ex, ey) * dyy(ex, ey) - dxy(ex, ey) ^ 2;
                curvature = (trace ^ 2) / determinant; 
                if (determinant < 0 || curvature >= ((r + 1) ^ 2) / r)
                    l_e(ex - 1, ey - 1) = 0;    % pay attention to '-1'
                end
            end
        end
        
        % gather all keypoints at current octave
        [x, y] = find(l_e);
        num_e = size(x);
        for p = 1 : num_e
            ex = x(p) + 1;
            ey = y(p) + 1;
            octave_kp = [octave_kp, [ex; ey]];
        end
        
    end
    % gather all keypoints at the target octave
    if (i == target_octave)
        X = octave_kp;
    end
end

S = target_octave;





