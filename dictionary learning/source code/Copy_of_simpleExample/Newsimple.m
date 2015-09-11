rng('default')

imgSize = 128;
B = 2; % Patch size
d = 1; % Patch skip (1=fully overlapping, B = nonoverlapping)
K = 128; % Dictionary size
pixelFraction = 0.05; % Portion of pixels to keep for inpainting
sigma = 0.05; % noise variance of added noise
numIter = 20;

fname = 'castle.png';
img0 = im2double(imread(fname));
img0 = imresize(img0, [imgSize,imgSize]); %image is actually 256x256
[N1,N2,N3] = size(img0);
noisy_img = img0 + sigma*randn(N1,N2, N3);


MaskSlice = binornd(1,pixelFraction, [N1,N2]);
sensingMask = repmat(MaskSlice, [1 1 N3]);
noisy_img = sensingMask.*noisy_img;


img1 = imresize(noisy_img, 0.5);
sensingMask1 = imresize(sensingMask,0.5);

img2 =imresize(noisy_img, 0.25);
sensingMask2 = imresize(sensingMask,0.25);

% ----------------------deal with the 0.25 -----------------

X = video2patches_fast(img2, B,B, d, d);
Phi = double(video2patches_fast(sensingMask2, B,B, d, d));

state = [];
for i = 1:20
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
end
D0 = state.D;
S0 = state.S;
Z0 = state.Z;
B0 = B;
state0 = state;


% expand dictionary
B = 2*B;
[D,S,Z] = vicEnlarge(D0,S0,B0,B);
state.D = D;
state.S = S;
state.Z = Z;

% ------------------deal with the 0.5----------------
X = video2patches_fast(img1, B,B, d, d);
Phi = double(video2patches_fast(sensingMask1, B, B, d, d));

for i = 1:10
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
end




% --------------------deal with 0.75 -------------
% B1 = 3 * B0;
% d = 3 * d0;
% state.D = vicResize(D0, B0, B1);
% 
% X = video2patches_fast(img10, B1,B1, d, d);
% Phi = double(video2patches_fast(sensingMask10, B1, B1, d, d));
% for i = 1:10
%     tic
%     state = BPFA_simple(X, K, Phi, state);
%     t = toc;
% end




% expand dictionary

B = 2*B;
[D,S,Z] = vicEnlarge(D0,S0,B0,B);
state.D = D;
state.S = S;
state.Z = Z;

% -----------------------deal with the 1.0----------------
X = video2patches_fast(noisy_img, B,B, d, d);
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

figure(1)
imagesc([img0, noisy_img, imgRecon]);
colormap('gray');
figure(2)
displayPatches(state.D);

T=toc;