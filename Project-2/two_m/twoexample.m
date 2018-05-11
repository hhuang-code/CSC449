
% load images using code in one_m
addpath(fullfile('..','one_m'));
oneload;

% a full example on an image
[S,C] = twosegment(15);
C = tworeduce(15,C,S);
B = twocut(S,C,61,0.5);
twoshow(15,C,S,B)
