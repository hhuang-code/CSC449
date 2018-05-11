function segments = tworeduce(index, segments, segmentimage)

% function segments = tworeduce (index,segments, segmentimage)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Wrapper for function to compute feature descriptors on the segments
%


global cards;

bins = 10;

im = cards{index};

for i=1:length(segments)

    segments(i).fv = histvec(im,segmentimage==i,bins);

end

