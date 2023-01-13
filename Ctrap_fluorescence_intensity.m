%% this script is to analyze the fluorescence intensity change over time of
% the C-trap scan data
mydir = '/Users/liuk3/Desktop/RM paper/datafigure/Scripts/example data files/'; % define the datapath
% first of all, I need to read the stacks
filename = 'Scan_movie_binding.tif';
ft = 5.23;  % define frame time
fname = [mydir filename];
Info = imfinfo(fname);
imageStack = [];
numframe = length(Info); % read the tiff stacks into matlab
for k = 1:numframe
    currentImage = imread(fname, k, 'Info', Info);
    imageStack(:,:,:,k) = currentImage;
end
%% now look at image
frames = 10; % define which frame it starts
currentImage = imageStack(:,:,:,frames);
channel = 1; % channel 1 is the LD650, and channel 2 is the SYTOX
I = currentImage(:,:,channel);
clim = [0 200]; % set up display range
close all;
imagesc(I,clim)
hold on;
roi = drawrectangle; % draw a rectangle roi
hold off
J = imcrop(I,roi.Position); % crop the image to the interest
rect = roi.Position; % define the crop region
imagesc(J,clim)
%% now generate the trace plot since it is moving I will do it for each 
close all;
% numframe = length(Info);
numframe = 29; % define the last frame to look
ph_mean=zeros(numframe-frames+1,1); % save the mean photon counts
channel = 1; % protein channel
for i = frames:numframe
    currentImage = imread(fname, i, 'Info', Info); % get the current frame image
    I = currentImage(:,:,channel); % get the protein channel
    J = imcrop(I,rect); % crop the image to the interest
    ph_mean(i-frames+1) = mean(mean(J)); % extract the mean photon counts
end
figure;
plot(ph_mean./ph_mean(1),'ro')
