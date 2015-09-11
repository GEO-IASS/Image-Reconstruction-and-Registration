rng('default')

imgSize = 128;
B = 12; % Patch size
d = 2; % Patch skip (1=fully overlapping, B = nonoverlapping)
K = 541; % Dictionary size
pixelFraction = 1; % Portion of pixels to keep for inpainting
sigma = 0.10; % noise variance of added noise
numIter = 40;

fname = 'house.png';%'/Users/d3y231/Downloads/gl.png';%
img0 = (im2double(imread(fname)));
%img0 = imresize(img0, [imgSize,imgSize]); %image is actually 256x256
[N1,N2,N3] = size(img0);

img = img0 + sigma*randn(N1,N2,N3);
sensingMask = binornd(1,pixelFraction, size(img));
img = sensingMask.*img;

X = video2patches_fast(img, B,B, d, d);
Phi = double(video2patches_fast(sensingMask, B,B, d, d));

fprintf('iter:\tpsnr\t\tsigmaEst\tS%%nonzero\tavgDUse\t\ttime\n');
state = [];
for i = 1:numIter
    tic
    state = BPFA_simple(X, K, Phi, state);
%     state = GHDL(X, K, Phi, state);
    t = toc;
    
    Xrecon = state.D*(state.S)';
    imgRecon = patches2video_fast(Xrecon, N1,N2,N3, d,d);
    
    fprintf('%4d:\t%f\t%f\t%f\t%f\t%f\n', i, psnr(img0,imgRecon), sqrt(1/state.geps), 100*nnz(state.S)/numel(state.S),full(mean(sum(state.S~=0,2))), t);
end

figure(1)
imagesc([img0, img, imgRecon]);
colormap('gray');
figure(2)
displayPatches(state.D);


