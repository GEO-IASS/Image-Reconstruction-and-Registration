tic;
tstart = tic;

% read the data and set the basic variables
load('C:\Users\haow889\Documents\MATLAB\data\gap160.mat');
[width,height,n] = size(x);
windowLength = 11;
iternum = n - windowLength + 1;
results = zeros(width,height,iternum);

% Generate sift features and store them in cells
siftinfo = cell(3,n);
for i = 1:n
    fname = sprintf('C:\\Users\\haow889\\Documents\\MATLAB\\data\\x_%d.png',i);
    [image, descriptors, locs] = sift(fname);
    s = struct('image',image,'des',descriptors,'locs',locs);
    siftinfo(:,i)=struct2cell(s);
end

% alignement block
DisRatio = 0.9;
p_value = 7; % p value ranges from 7 to 20 
for i = 1:iternum
    % use window to slide the whole image sets and for each slide step, fix the
    % middle one as the reference match the former half and latter half
    alignedImages = zeros (width, height, windowLength);
    tempMiddle = i + (windowLength-1)/2;
    siftinfo1 = siftinfo(:,tempMiddle);
    alignedImages(:,:,round(windowLength/2)) = cell2mat(siftinfo1(1,1));
    des1 = cell2mat(siftinfo(2,tempMiddle));
    loc1 = cell2mat(siftinfo(3,tempMiddle));
    
    % match the former half
    for j = i : tempMiddle-1
        siftinfo2 = siftinfo(:,j);
        alignedImages(:,:,j-i+1)=alignment(siftinfo1,siftinfo2,DisRatio);
    end
    
    % match the latter half
    for j = tempMiddle+1 : i+windowLength-1
        siftinfo2 = siftinfo(:,j);
        alignedImages(:,:,j-i+1)=alignment(siftinfo1,siftinfo2,DisRatio);
    end
        
    
    
    % aggregation block
    result = aggregation(alignedImages,p_value);
    results(:,:,i) = (result-min(result(:)))./(max(result(:))); % normalize
    resultName = sprintf('alignagre%d_%d.png',i,i+windowLength-1);
    imwrite(results(:,:,i),strcat('C:\Users\haow889\Documents\MATLAB\',...
        'figure\alignagre_images\',resultName),'png');
end
telapsed = toc(tstart);

