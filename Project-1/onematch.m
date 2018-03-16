function M = onematch(F1,F2,k)

% function M = onereduce(F1,F2,k)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Wrapper for function to matching extracted feature vectors from a pair
%  of images
%
%  F1 is the feature matrix (rows -dimensions and cols number of points)
%       from image 1
%  F2 feature matrix from image 2
%  k is the number of matches to take (optional)
%
%  M is a k x 2 matrix where k is the number of matches and the first col
%    is the index of the match in F1 and the second col in F2

if nargin==2
    k=10;
end

n1 = size(F1,2);
n2 = size(F2,2);

D = ones(n1,1);
I = ones(n1,1);

for i=1:n1
   % I is the matched index of F2; D is distance
   [I(i),D(i)] = onematchsingle(F1(:,i),F2);
end

% now, rank and take just the top $k=5$
[Ds,Di] = sort(D,'ascend'); % Ds is distance; Di is the matched index of F1

M=zeros(k,2);
for i=1:k
    M(i,1) = Di(i);
    M(i,2) = I(Di(i));
end
