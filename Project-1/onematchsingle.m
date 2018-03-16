function [m,d] = onematchsingle(f,F)

% function [m,d] = onematchsingle(f,F)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Wrapper for function to matching an feature vector to a feature matrix
%
%  f is the vector 
%  F is the matrix, each column is a feature vector
%
% m is the matched index
% d is the distance for the match


%%%%%%%%%%% add your code below

dist = sqrt(sum(bsxfun(@minus, F, f) .^ 2, 1));
[d, m] = min(dist);  % pay attention to the order of m and d

%%%%%%%%%% add your code above
