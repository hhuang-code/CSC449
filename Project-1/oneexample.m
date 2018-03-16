% an example of a full run through two images toward display

clear

if ~(exist('cards','var')==1)
    oneload2
end
if ~(exist('index1','var')==1)
    index1 = 1;
end
if ~(exist('index2','var')==1)
    index2 = 2;
end

[F1,X1] = onereduce(index1);
[F2,X2] = onereduce(index2);
M = onematch(F1,F2);
oneshow(index1,X1,index2,X2,M);
saveas(gcf,'lastexample.png','png');