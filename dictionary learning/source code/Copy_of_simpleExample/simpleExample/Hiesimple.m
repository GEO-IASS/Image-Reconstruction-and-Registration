tic;

rng('default');

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

% use the same sensingMask, imresize the noisy_img but smooth off the
% black holes of the original image

% img1 = imresize(noisy_img, 0.5);
% sensingMask1 = imresize(sensingMask, 0.5);
% 
% img2 = imresize(noisy_img0, 0.25);
% sensingMask2 = imresize(sensingMask, 0.25);
img1 = imresize(img0, 0.5);
[N11, N12, N13] = size(img1);
img1 = img1 + sigma*randn(N11,N12, N13);
MaskSlice1 = binornd(1,pixelFraction, [N11,N12]);
sensingMask1 = repmat(MaskSlice1, [1 1 N13]);
% sensingMask1 = imresize(sensingMask, 0.5);
img1 = sensingMask1.*img1;

img2 = imresize(img0, 0.25);
[N21, N22, N23] = size(img2);
img2 = img2 + sigma*randn(N21,N22, N23);
MaskSlice2 = binornd(1,pixelFraction, [N21,N22]);
sensingMask2 = repmat(MaskSlice2, [1 1 N23]);
% sensingMask2 = imresize(sensingMask, 0.25);
img2 = sensingMask2.*img2;
%{
% generate new sensingMask,imresize the original image and 
% re-sample the resized images

% img10 = imresize(img0, 0.75);
% [N101, N102, N103] = size(img10);
% img10 = img10 + sigma*randn(N101,N102, N103);
% MaskSlice10 = binornd(1,pixelFraction, [N101,N102]);
% sensingMask10 = repmat(MaskSlice10, [1 1 N103]);
% img10 = sensingMask10.*img10;



img1 = imresize(img0, 0.5);
[N11, N12, N13] = size(img1);
img1 = img1 + sigma*randn(N11,N12, N13);
MaskSlice1 = binornd(1,pixelFraction, [N11,N12]);
sensingMask1 = repmat(MaskSlice1, [1 1 N13]);
img1 = sensingMask1.*img1;


% img21 = imresize(img0, 0.375);
% [N211, N212, N213] = size(img21);
% img21 = img21 + sigma*randn(N211,N212, N213);
% MaskSlice21 = binornd(1,pixelFraction, [N211,N212]);
% sensingMask21 = repmat(MaskSlice21, [1 1 N213]);
% img21 = sensingMask21.*img21;



img2 = imresize(img0, 0.25);
[N21, N22, N23] = size(img2);
img2 = img2 + sigma*randn(N21,N22, N23);
MaskSlice2 = binornd(1,pixelFraction, [N21,N22]);
sensingMask2 = repmat(MaskSlice2, [1 1 N23]);
img2 = sensingMask2.*img2;

% img3 = imresize(img0, 0.125);
% [N31, N32, N33] = size(img3);
% img3 = img3 + sigma*randn(N31,N32, N33);
% MaskSlice3 = binornd(1,pixelFraction, [N31,N32]);
% sensingMask3 = repmat(MaskSlice3, [1 1 N33]);
% img3 = sensingMask3.*img3;


% img4 = imresize(img0, 0.0625);
% [N41, N42, N43] = size(img4);
% img4 = img4 + sigma*randn(N41,N42, N43);
% MaskSlice4 = binornd(1,pixelFraction, [N41,N42]);
% sensingMask4 = repmat(MaskSlice4, [1 1 N43]);
% img4 = sensingMask4.*img4;
%}
% imresize the original image and create the corresponding sensingMask

%-------------------------deal with 0.0625--------------
%{
X = video2patches_fast(img4, B,B, d, d);
Phi = double(video2patches_fast(sensingMask4, B,B, d, d));
state = [];
for i = 1:30
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
end
D0 = state.D;
B0 = B;
d0 = d;

% expand dictionary
% Note it is not exactly ABCD -> AAAABBBBCCCCDDDD
B = 2*B;
d = 2*d;
D = zeros(B^2*N3, K);
for j = 1:K
    element = D0(:,j);
%     temp = imresize(element,4);
%     D(:,j) = temp(:,1);
    temp = myresize(element,B/B0,B0,d0);
    D(:,j) = temp;
end
state.D = D;

%}

% ----------------------deal with the 0.125 -----------------
%{
X = video2patches_fast(img3, B,B, d, d);
Phi = double(video2patches_fast(sensingMask3, B,B, d, d));

% state = [];
for i = 1:20
    tic
    state = BPFA_simple(X, K, Phi, state);
    t = toc;
end
% D0 = state.D;
% B0 = B;
% d0 = d;

% expand dictionary
% Note it is not exactly ABCD -> AAAABBBBCCCCDDDD
B = 2*B;
d = 2*d;
D = zeros(B^2*N3, K);
for j = 1:K
    element = D0(:,j);
%     temp = imresize(element,4);
%     D(:,j) = temp(:,1);
    temp = myresize(element,B/B0,B0,d0);
    D(:,j) = temp;
end
state.D = D;
%}

% ----------------------deal with the 0.25 -----------------

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
D0 = state.D;
B0 = B;
d0 = d;


% expand dictionary
B1 = 2 * B0;
d1 = 2 * d0;
D = vicResize(D0,B0,B1);
% D = zeros(B1^2*N3, K);
% for j = 1:K
%     element = D0(:,j);
%     cube0 = patches2video_fast(element,B0,B0,3,d0,d0);
%     temp = imresize(cube0,B1/B0);
%     D(:,j) = video2patches_fast(temp,B1,B1,d1,d1);
% %     temp = myresize(element,B/B0,B0,d0);
% %     D(:,j) = temp;
% end
state.D = D;

% B1 = 2 * B0;
% d = 2 * d0;
% state.D = vicResize(D0, B0, B1);



% ------------------deal with the 0.5----------------
X = video2patches_fast(img1, B1,B1, d1, d1);
Phi = double(video2patches_fast(sensingMask1, B1, B1, d1, d1));
% resize the S and Z

% [P,N] = size(X);
% state.S = sprandn(N,K,0.2);
% state.Z = full(state.S~=0);
% state.Pi = mean(state.Z);


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

B1 = 2*B1;
d1 = 2*d1;
D = vicResize(D0,B0,B1);
% D = zeros(B1^2*N3, K);
% for j = 1:K
%     element = D0(:,j);
%     cube0 = patches2video_fast(element,B0,B0,3,d0,d0);
%     temp = imresize(cube0, B1/B0);
%     D(:,j) = video2patches_fast(temp,B1,B1,d1,d1);
% %     temp = myresize(element,B/B0,B0,d0);
% %     D(:,j) = temp;
% end
state.D = D;
% B1 = 4 * B0;
% state.D = vicResize(D0, B0, B1);

% -----------------------deal with the 1.0----------------
X = video2patches_fast(noisy_img, B1, B1, d1, d1);
Phi = double(video2patches_fast(sensingMask, B1, B1, d1, d1));

% resize the S and Z
% [P,N] = size(X);
% state.S = sprandn(N,K,0.2);
% state.Z = full(state.S~=0);
% state.Pi = mean(state.Z);


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