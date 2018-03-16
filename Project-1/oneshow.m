function oneshow(index1,X1,index2,X2,M)

% function oneshow(index1,X1,index2,X2,M)
%
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
%  M is a k x 2 matrix where k is the number of matches and the first col
%    is the index of the match in F1 and the second col in F2

global cards;

% just display detected points on the image
if nargin==2
    
    figure; 
    imagesc(cards{index1});
    hold on;
    
    plot(X1(1,:),X1(2,:),'r+');
    
    title('oneshow: detected points');
    
    return;
end

% display the points in two images side by side
if nargin==4
    figure;
    subplot(1,2,1);
    imagesc(cards{index1});
    hold on;
    plot(X1(1,:),X1(2,:),'r+');
    
    subplot(1,2,2);
    imagesc(cards{index2});
    hold on;
    plot(X2(1,:),X2(2,:),'r+');
    
    return;
end

% display the matches too
if nargin==5
    
    figure;
    im = [cards{index1},cards{index2}];
    imagesc(im);
    hold on;
    [r,c,b] = size(cards{index1});
    
    plot(X1(1,:),X1(2,:),'r+');
    plot(X2(1,:)+c,X2(2,:),'r+');
    
    for i=1:size(M,1)
        line([X1(1,M(i,1));X2(1,M(i,2))+c],  ...
             [X1(2,M(i,1));X2(2,M(i,2))],'color',rand(1,3),'linewidth',2);
    end
    
    return;
end

fprintf('oneshow does not understand the arguments. see help oneshow\n');


