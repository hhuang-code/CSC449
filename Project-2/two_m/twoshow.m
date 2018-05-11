function twoshow(index,centers,segmentimage,mask)

% function twoshow(index,centers,segmentimage)
%
%     Jason Corso
%     Chenliang Xu
%
% different behaviors with 2, 3 or 4 parameters.  


global cards;

if nargin > 3
    if (~isempty(mask))
        figure;
        I = cards{index};
        for i=1:3
            Ii = I(:,:,i);
            Ii(mask~=1) = 0;
            I(:,:,i) = Ii;
        end
        
        imagesc(I);
    end
end

if nargin > 2
    if ~isempty(segmentimage)
        
        dy = conv2(segmentimage,fspecial('sobel'),'same');
        dx = conv2(segmentimage,fspecial('sobel')','same');
        M = (dx.*dx+dy.*dy)>0;
        I = cards{index};
        for i=1:3
            Ii = I(:,:,i);
            Ii(M) = 0;
            I(:,:,i) = Ii;
        end
        
        figure;
        subplot(1,2,1);
        imagesc(I);
        subplot(1,2,2);
        imagesc(segmentimage);
        
        return
    end
end

figure;
imagesc(cards{index});
hold on;

if ~isempty(centers)
    for i=1:length(centers)
        plot(centers(i).x,centers(i).y,'g+');
    end
end


