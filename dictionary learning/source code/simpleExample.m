tic;
rng('default')

imgSize = 128;
B = 8; % Patch size
d = 1; % Patch skip (1=fully overlapping, B = nonoverlapping)
K = 128; % Dictionary size
pixelFraction = 0.05; % Portion of pixels to keep for inpainting
sigma = 0.05; % noise variance of added noise
numIter = 5;

fname = 'castle.png';
img0 = im2double(imread(fname));
img0 = imresize(img0, [imgSize,imgSize]); %image is actually 256x256
[N1,N2,N3] = size(img0);

img = img0 + sigma*randn(N1,N2, N3);
sensingMask = zeros(N1,N2,N3);
MaskSlice = binornd(1,pixelFraction, [N1,N2]);
for i = 1:N3
    sensingMask(:,:,i) = MaskSlice;
end
img = sensingMask.*img;



X = video2patches_fast(img, B,B, d, d);
Phi = double(video2patches_fast(sensingMask, B,B, d, d));

psnrBefore = psnr(img0, img);
fprintf('iter:\tpsnr\t\tsigmaEst\tAvgDUse\t\ttime\n');
state = [];
for i = 1:numIter
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
    
    Xrecon = state.D*(state.Z.*state.S)';
    imgRecon = patches2video_fast(Xrecon, N1,N2,N3, d,d);
    
    fprintf('%4d:\t%f\t%f\t%f\t%f\n', i, psnr(img0,imgRecon), sqrt(1/state.geps), mean(sum(state.Z,2)), t);
end
    
psnrAfter = psnr(img0, imgRecon);

figure(1)
imagesc([img0, img, imgRecon]);
colormap('gray');
figure(2)
displayPatches(state.D);

T=toc;
