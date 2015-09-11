rng('default')
data = load('C:\Users\haow889\Documents\MATLAB\data\specImage.mat');


imgSize = 48;
B = 2; % Patch size
d = 1; % Patch skip (1=fully overlapping, B = nonoverlapping)
K = 12; % Dictionary size
pixelFraction = 0.5; % Portion of pixels to keep for inpainting
sigma = 0.05; % noise variance of added noise
numIter = 20;


img0 = im2double(data.img);
img0 = imresize(img0, [imgSize,imgSize]); %image is actually 256x256
[N1,N2,N3] = size(img0);

% normalize each slice 
for i = 1:N3
    img0(:,:,i) = (img0(:,:,i) - min(min(min(img0(:,:,i)))))./max(max(max(img0(:,:,i))));
end

noisy_img = img0 + sigma*randn(N1,N2, N3);
sensingMask = zeros(N1,N2,N3);
MaskSlice = binornd(1,pixelFraction, [N1,N2]);
for i = 1:N3
    sensingMask(:,:,i) = MaskSlice;
end
noisy_img = sensingMask.*noisy_img;


img1 = imresize(noisy_img, 0.5);
sensingMask1 = sensingMask;
sensingMask1 = imresize(sensingMask1, 0.5);

img2 =imresize(noisy_img, 0.25);
sensingMask2 = sensingMask;
sensingMask2 = imresize(sensingMask2, 0.25);


% deal with the 0.25 

X = video2patches_fast(img2, B,B, d, d);
Phi = double(video2patches_fast(sensingMask2, B,B, d, d));

% psnrBefore = psnr(img0, img);
% fprintf('iter:\tpsnr\t\tsigmaEst\tAvgDUse\t\ttime\n');
state = [];
for i = 1:20
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
end
B0 = B;
d0 = d;
% expand dictionary
% Note it is not exactly ABCD -> AAAABBBBCCCCDDDD
B = 2*B;
d = 2*d;
D = zeros(B^2*N3, K);
for j = 1:K
    element = state.D(:,j);
    cube0 = reshap(element,B0,B0,1,1);
    temp = imresize(element, 4);
    D(:,j) = temp(:,1);
    
end
state.D = D;



% deal with the 0.5
X = video2patches_fast(img1, B,B, d, d);
Phi = double(video2patches_fast(sensingMask1, B, B, d, d));
for i = 1:10
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
end

% expand dictionary
% Note it is not exactly ABCD -> AAAABBBBCCCCDDDD
B = 2*B;
d = 2*d;
D = zeros(B^2*N3, K);
for j = 1:K
    element = state.D(:,j);
    temp = imresize(element, 4);
    D(:,j) = temp(:,1);
end
state.D = D;

% deal with the 1.0
X = video2patches_fast(noisy_img, B, B, d, d);
Phi = double(video2patches_fast(sensingMask, B, B, d, d));
for i = 1:5
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
    
    Xrecon = state.D*(state.Z.*state.S)';
    imgRecon = patches2video_fast(Xrecon, N1,N2,N3, d,d);
    fprintf('%4d:\t%f\t%f\t%f\t%f\n', i, psnr(img0,imgRecon), sqrt(1/state.geps), mean(sum(state.Z,2)), t);
end

   
psnrAfter = psnr(img0, imgRecon);




