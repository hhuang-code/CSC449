function [l,Copt] = slic(image,k)

% function [S,Copt] = slic(image,numsegs)
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% Implementation of Achanta et al. "SLIC Superpixels Compared 
%  to State-of-the-art Superpixel Methods" PAMI 2011.
% Notation is directly from this paper when possible.
%
% k is the number of segments
% l is an index image with a unique id associated with each segments
% Copt is the structure array describing the segment information
%    it is an optional output


[r,c,b] = size(image);
N = r*c;
S = ceil(sqrt(N/k));

fprintf('r=%d,c=%d,N=r*c=%d\n',r,c,N);
fprintf('S=%d\n',S);

% convert the image to Lab space
I = RGB2Lab(image);
IL = I(:,:,1);
Ia = I(:,:,2);
Ib = I(:,:,3);
dy = conv2(IL,fspecial('sobel'),'same');
dx = conv2(IL,fspecial('sobel')','same');

% set up how we will store the cluster centers
C = struct( 'index',{}...         
           ,'l'    ,{}...
           ,'a'    ,{}...
           ,'b'    ,{}...
           ,'x'    ,{}...
           ,'y'    ,{}...
           ,'x_sub',{}... % subpixel x
           ,'y_sub',{}... % subpixel y
           ,'fv'   ,{}... % a feature vector that will get used later
          );
            
% initialize the seed points
i=1;
for y=S/2:S:r
    for x=S/2:S:c
        C(i).index = i;
        C(i).x = x;
        C(i).y = y;
        C(i).l = IL(y,x);
        C(i).a = Ia(y,x);
        C(i).b = Ib(y,x);
        i = i + 1;
    end
end
clear i x y;

% assert length(centers) = k
k = length(C);

% recenter seed points
M = dx.*dx+dy.*dy;
for i = 1:k
    minr = max(1,C(i).y-1);
    minc = max(1,C(i).x-1);
    maxr = min(r,C(i).y+1);
    maxc = min(c,C(i).x+1);
    g = M(minr:maxr,minc:maxc);
    
    [~,ix] = min(g(:));
    [y,x] = ind2sub([3 3],ix);
    y = C(i).y+y-2;
    x = C(i).x+x-2;
    C(i).x = x;
    C(i).y = y;
    C(i).l = IL(y,x);
    C(i).a = Ia(y,x);
    C(i).b = Ib(y,x);
end

% initialize other variables
% label l
l = ones(r,c)*-1;   % all -1 matrix
% distance
d = inf(r,c);
% residual error
E = 1e20;
% m -- scalar weight on the distance function
m = 15;
m2overS2 = m*m/(S*S);   % m^2 / S^2

% main loop of algorithm
for iteration=1:12   % max 12 iterations
    
    %%%%%  fill in the main body of the SLIC algorithm here
    for i = 1 : k   % for each cluster
        for y = max(1, C(i).y - S) : min(r, C(i).y + S - 1)
            for x = max(1, C(i).x - S) : min(c, C(i).x + S - 1)
                % for each pixel in 2S * 2S regin around this cluster
                % fprintf('i=%d, y=%d, x=%d, C(i).x=%d, C(i).y=%d\n',i,y,x,C(i).x,C(i).y);
                dc = (IL(y, x) - IL(C(i).y, C(i).x))^2 + (Ia(y, x) - Ia(C(i).y, C(i).x))^2 + (Ib(y, x) - Ib(C(i).y, C(i).x))^2;
                ds = (y - C(i).y)^2 + (x - C(i).x)^2;
                dist = sqrt(dc + ds * m2overS2);
                if dist < d(y, x)
                    d(y, x) = dist;
                    l(y, x) = i;
                end
            end
        end
    end
    
    % at each iteration, compute the current residual error (sum of change in x,y for all segments)
    %  use a variable named err for that.  
    
    NC = C;
    % reset NC (new clusters)
    for i = 1 : k
        NC(i).x = 0;
        NC(i).y = 0;
        NC(i).l = 0;
        NC(i).a = 0;
        NC(i).b = 0;
    end
    
    % compute new value of NC
    for y = 1 : r
        for x = 1 : c
            NC(l(y, x)).y = NC(l(y, x)).y + y;
            NC(l(y, x)).x = NC(l(y, x)).x + x;
            NC(l(y, x)).l = NC(l(y, x)).l + IL(y, x);
            NC(l(y, x)).a = NC(l(y, x)).a + Ia(y, x);
            NC(l(y, x)).b = NC(l(y, x)).b + Ib(y, x);
        end
    end
    
    % compute new center
    for i = 1 : k
        tot = sum(sum(l == i));
        NC(i).y = round(NC(i).y / tot);
        NC(i).x = round(NC(i).x / tot);
        NC(i).l = NC(i).l / tot;
        NC(i).a = NC(i).a / tot;
        NC(i).b = NC(i).b / tot;
    end
    
    % compute residual error
    err = 0;
    for i = 1 : k
        err = err + sqrt((NC(i).y - C(i).y)^2 + (NC(i).x - C(i).x)^2 + (NC(i).l - C(i).l)^2 + (NC(i).a - C(i).a)^2 + (NC(i).b - C(i).b)^2);
    end
    
    C = NC;
    
    %%%%% do not change below this line

    % check residual
    if (E-err) < 1
       break;
    else
        E = err;
    end
    
end


% finish
if nargout==2
    Copt = C;
end
