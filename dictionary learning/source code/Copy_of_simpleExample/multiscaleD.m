rng('default')

imgSize = 128;
B0 = 12; % Patch size
d = 2; % Patch skip (1=fully overlapping, B = nonoverlapping)

pixelFraction = 1; % Portion of pixels to keep for inpainting
sigma = 0.10; % noise variance of added noise

fname = 'castle.png';%'/Users/d3y231/Downloads/gl.png';%
img0 = rgb2gray(im2double(imread(fname)));
%img0 = imresize(img0, [imgSize,imgSize]); %image is actually 256x256
[N1,N2,N3] = size(img0);

img = img0 + sigma*randn(N1,N2,N3);
masks{1} = binornd(1,pixelFraction, size(img));
imgs{1} = masks{1}.*img;

imgs{2} = imresize(imgs{1},0.5); masks{2} = imgs{2} ~= 0;
imgs{3} = imresize(imgs{2},0.5); masks{3} = imgs{3} ~= 0;

B = B0;
X{1} = video2patches_fast(imgs{1}, B,B, d, d);
Phi{1} = double(video2patches_fast(masks{1}, B,B, d, d));
B = B/2;
X{2} = video2patches_fast(imgs{2}, B,B, d, d);
Phi{2} = double(video2patches_fast(masks{2}, B,B, d, d));
B = B/2;
X{3} = video2patches_fast(imgs{3}, B,B, d, d);
Phi{3} = double(video2patches_fast(masks{3}, B,B, d, d));
% 
szs =    [1,2,3  ,2,1];
scales = [1,0.5,0.5,2,2,1];
numIter = [2,5,10,5,40];

% szs =    [2,3  ,2,1];
% scales = [1,0.5,2,2,1];
% numIter = [2,2,2,5,5];

K = 5; % Dictionary size
[P,N] = size(X{szs(1)});
state.D = randn(P,K);
state.S = sprandn(N,K,0.2);
state.Z = full(state.S~=0);
state.Pi = mean(state.Z);
state.gs = ones(1,K);
state.geps = 1;

for i = 1:length(szs)
    [P,N] = size(X{szs(i)});
    [~,K] = size(state.D);
%     state.Z = logical(binornd(1,repmat(state.Pi,N,1)));
%     state.S = (pinv(state.D)*X{szs(i)} .* state.Z')';
    state.S = sprandn(N,K,0.2);
    state.Z = full(state.S~=0);
    state.Pi = mean(state.Z);
    state.gs = ones(1,K);
    for j = 1:numIter(i)
        tic
        state = BPFA_simple(X{szs(i)}, K, Phi{szs(i)}, state);
        t = toc;

        Xrecon = state.D*(state.S)';

        fprintf('%4d-%4d:\t%f\t%f\t%f\n', i,j, mse(Xrecon, X{szs(i)}), sqrt(1/state.geps), t);
    end
    figure(2+szs(i)-1)
    displayPatches(state.D);drawnow;
    
    state = resizeD(state, scales(i+1), 1, 1);
end

imgRecon = patches2video_fast(Xrecon, N1,N2,N3, d,d);

figure(1)
imagesc([img0, imgs{1}, imgRecon]);
colormap('gray');



