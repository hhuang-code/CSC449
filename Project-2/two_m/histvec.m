function v = histvec(image,mask,b)

% function v = histvec(image,mask,b)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
%  For each channel in the image, compute a b-bin histogram (uniformly space
%  bins in the range 0:1) of the pixels in image where the mask is true. 
%  Then, concatenate the vectors together into one column vector (first
%  channel at top).
%
%  mask is a matrix of booleans the same size as image.
% 
%  normalize the histogram of each channel so that it sums to 1.
%
%  You CAN use the hist function.  (since you have already worked on 
%    implementing hist from assignment 1)

chan = size(image,3);

c = 1/b;      % bin offset
x = c/2:c:1;  % bin centers


%%%%% implement below this line into a 3*b by 1 vector  v
%%  3*b because we have a color image and you have a separate 
%%  histogram per color channel

v = [];

for i = 1 : chan    % for each channel   
    single = image(:,:,i);  % split channels
    ch = single(mask);  % apply mask
    h = hist(ch, x);    % compute histogram
    v = cat(1, v, h' ./ sum(h));
end

v = v';  % from row vector to column vector

