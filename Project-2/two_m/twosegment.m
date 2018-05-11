function [S,C] = twosegment(index, numsegs)

% function [S,C] = twosegment(index,numsegs)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Wrapper for function to compute segments/superpixels on an image
%
%  numsegs is optional and is the number of segments desired
%
% S is an index image with a unique id associated with each segments
% C is a struct-array describing the segments themselves

global cards;

if nargin==1
    numsegs = 128;
end

im = cards{index};

[S,C] = slic(im,numsegs);


