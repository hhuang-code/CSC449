function Bmap = segNeighbors(svMap)

% function segNeighbors(svMap)
%
%     Chenliang Xu
%
% Return the adjacency map of segmentation

segmentList = unique(svMap);
segmentNum = length(segmentList);

Bmap = zeros(segmentNum);
for i=1:length(segmentList)
    svmap_i = (svMap==segmentList(i));
    svmap_i_b = imdilate(svmap_i, ones(3,3))&(~svmap_i);
    neighbors = unique(svMap(svmap_i_b));
    Bmap(segmentList(i), neighbors) = 1;
end
Bmap = sparse(Bmap); % optional
