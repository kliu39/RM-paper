%% this is a matlab script to analyze the contour length of DNA molecules from 
% AFM images. 
datapath =  '/Users/liuk3/Desktop/RM paper/datafigure/Scripts/example data files/';
filename = 'DNA_RM_0000.ibw';
D = IBWread([datapath filename]); % read the ibw file, make sure the function IBWread.m is in path.
layer = 1;
I  = D.y(:,:,layer); % read as image file
angle = 90;
J = imrotate(I,angle); % the default rotates 90 degrees
A = size(J);
B = reshape(J,[A(1)*A(2),1]);    % convert this to 1D array
perline = 0.98; % what is the percentage to filter the image, adjust accordingly
B = sort(B);
level = B(round(length(B)*perline)); % find the actual height value to filter
BW = imbinarize(J,level); % apply filter to the new image to binarize
bskel = bwskel(BW,'MinBranchLength',0); % build the skeleton
L = bwlabel(bskel); % label the image, L is the same size of the binary image with labeled number at the positions. 
imshowpair(I,BW,'montage') % show the original and the binarized images
%% now trace and measure the length
close all; 
bskel = bwskel(BW,'MinBranchLength',0); % build the skeleton
EP = find(bwmorph(bskel, 'endpoints')); % find the endpoints of the branches; the order of EP is the same as epX and Y. 
EP(:,2) = L(EP); % label the end points
BP = find(bwmorph(bskel, 'branchpoints')); % find the branchpoints of the branches
BP(:,2) = L(BP); % label the branch points
[epX, epY] = find( bwmorph(bskel, 'endpoints') ); % get the X and Y coordinates directly.
[bpX, bpY] = find( bwmorph(bskel, 'branchpoints') ); % get the X and Y coordinates directly.
[~,ia] = unique(EP(:,2)); % get the unique index in the EP list.
% now I need to label the data how many ends it has.
N_ep = zeros(length(ia),1); % number of ep
EtEd = zeros(length(ia),1); % end to end distance

figure % plot the labeled results
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on;
trace = cell(length(ia),1);
trace_L = L(EP(ia));
for i = 1 : length(ia)
    trace{i} = bwtraceboundary(bskel, [epX(ia(i)), epY(ia(i))],'E');    
    plot(trace{i}(:,2),trace{i}(:,1),'g','LineWidth',2)
    text(epY(ia(i)),epX(ia(i)), sprintf('%i',i),'fontsize',FS)
end

figure;
imagesc(J)
%% now calculate the length based on label
t_len = zeros(length(trace),1); % store the total trace length
S = struct(); % store segment length with branchpoints
ed_len = struct(); % store end to end distance of endpoints.
edge = 30; % define pixel numbers from edge to exclude
ifilt = size(J,1)-edge > epX(ia) & epX(ia) > edge & epY(ia) > edge & size(J,2)-edge > epY(ia); % filter out the edge
% idx_f = ia(ifilt);
idx_f = find(ifilt==1); % only look for traces away from the edge 
for i = 1 : length(idx_f)
    ln = EP(ia(idx_f(i)),2); % label number
    ep = find(EP(:,2)==ln); % find the corresponding end points
    bp = find(BP(:,2)==ln); % find the corresponding branch points
    if length(ep)>1
        idx = 1:(length(trace{idx_f(i)})-1)/2; % trace to the half point will be the full length
        diff_1 = diff(trace{idx_f(i)}(idx,1))*D.dx(1)*10^9;  % the difference matrix of dimension 1
        diff_2 = diff(trace{idx_f(i)}(idx,2))*D.dx(2)*10^9; % the difference matrix of dimension 2
        t_len(idx_f(i)) = sum(sqrt(diff_1.^2+diff_2.^2)); % total contour length
    end
    if ~isempty(bp) % if there are branch points
        S(i).label = ln;
        for k = 1 : length(bp) % for each branch points, calculate it is distance to the ends.
            xbp = bpX(bp(k),1); % break points X coordinates
            ybp = bpY(bp(k),1); % break points Y coordinates
            bp_idx = find(trace{idx_f(i)}(:,1)==xbp & trace{idx_f(i)}(:,2)==ybp); % find the segment
            for m = 1 : length(ep) % end point
                xep = epX(ep(m),1);
                yep = epY(ep(m),1);
                ep_idx = find(trace{idx_f(i)}(:,1)==xep & trace{idx_f(i)}(:,2)==yep);
                idx_d = sort([ep_idx(1),bp_idx(1)]);
                idx = idx_d(1):idx_d(2);
                diff_1 = diff(trace{idx_f(i)}(idx,1))*D.dx(1)*10^9;
                diff_2 = diff(trace{idx_f(i)}(idx,2))*D.dx(2)*10^9;
                S(i).s_len(k,m) = sum(sqrt(diff_1.^2+diff_2.^2)); % calculate the segment length for each DNA
            end
        end
    end
end