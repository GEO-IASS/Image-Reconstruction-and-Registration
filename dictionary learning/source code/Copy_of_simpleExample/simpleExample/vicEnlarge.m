function [D1,S1,Z1] = vicEnlarge(D,S,B,B1)
    [M, N] = size(D);
    slicenum = M / B^2;
    D1 = zeros(B1^2*slicenum, N);
    for i = 1:N
        vector = D(:,i);
        cube = patches2video_fast(vector,B,B,slicenum,1,1);
        cube1 = imresize(cube,B1/B);
        D1(:,i) = video2patches_fast(cube1,B1,B1,1,1);
    end
    clear M N;
    [M,N] = size(S);
    colD = sqrt(M);
    index = ones(colD,colD);
   
    imsize1 = B1/B *(colD + B - 1);
    colD1 = imsize1 - B1 + 1;
    index1(1:B1/B:colD1,1:B1/B:colD1) = index;
    index1_1D = video2patches_fast(index1,1,1,1,1);
    M1 = colD1^2;
    S1 = zeros(M1,N);
    count = 1;
    for i = 1:M1
        if index1_1D(1,i) == 1
            S1(i,:) = S(count,:);
            count = count + 1;
        else
            S1(i,:) = sprandn(1,N,0.2);
        end
    end
    [i,j,s] = find(S1);
    S1 = sparse(i,j,s,M1,N);
    Z1 = full(S1~=0);
 
     
end