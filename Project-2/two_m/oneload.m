% oneload.m
% 
%     EECS Foundation of Computer Vision;
%     Jason Corso
%
% script to load all of the images in the dataset
%
global cards;


D = dir(fullfile('cards','*png'));

fprintf('loading %02d',length(D));

cards = cell(20,1);

for i = 1:length(D)
    fprintf('\b\b%02d',i);
    cards(i) = {double(imread(fullfile('cards',D(i).name)))/255.0};
end
fprintf('\n');

clear D;
