function [F,Xout] = onereduce(index)

% function [F,Xout] = onereduce(index)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Wrapper for function to detect and describe interest points in image
%
%  index is the index for the image
%
%  output
%  F is a kxn matrix of reduced features
%         n is the number of features
%         k is the dimensionality of each feature
%  Xout (optional) are the feature point locations
%


global cards;

X = onedetect(index);
n = size(X,2);

F = cell(1,n);

for i = 1:n
   F(1,i) = {hog(cards{index},X(1,i),X(2,i),16)}; 
end

F = cell2mat(F);

if nargout == 2
    Xout = X;
end
