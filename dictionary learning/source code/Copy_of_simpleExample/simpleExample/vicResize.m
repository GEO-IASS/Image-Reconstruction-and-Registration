function D1 = vicResize(D,B,B1)
    [M, N] = size(D);
    slicenum = M / B^2;
    D1 = zeros(B1^2*slicenum, N);
    for i = 1:N
        vector = D(:,i);
        cube = patches2video_fast(vector,B,B,slicenum,1,1);
        cube1 = zeros(B1,B1,slicenum);
        if mod(B1,2) == 0
            num = B1 / B;
            cube1(1:B1/2, 1:B1/2, :)= repmat(cube(1,1,:),[num num 1]);
            cube1(1:B1/2, B1/2+1:B1, :) = repmat(cube(1,2,:),[num num 1]);
            cube1(B1/2+1:B1, 1:B1/2, :) = repmat(cube(2,1,:),[num num 1]);
            cube1(B1/2+1:B1, B1/2+1:B1,:) = repmat(cube1(2,2,:),[num num 1]);
            D1(:,i) = video2patches_fast(cube1,B1,B1,1,1);
        else
            B1 = B1 - 1;
            num = B1 / B;
            cube1(1:B1/2, 1:B1/2, :)= repmat(cube(1,1,:),[num num 1]);
            cube1(1:B1/2, B1/2+1:B1, :) = repmat(cube(1,2,:),[num num 1]);
            cube1(B1/2+1:B1, 1:B1/2, :) = repmat(cube(2,1,:),[num num 1]);
            cube1(B1/2+1:B1, B1/2+1:B1,:) = repmat(cube1(2,2,:),[num num 1]);
            for j = 1:B1
                cube1(j,B1+1,:) = cube1(j,B1,:);
                cube1(B1+1,j,:) = cube1(B1,j,:);
            end
            cube1(B1 + 1, B1 + 1,:) = cube1(B1, B1, :);
            B1 = B1 + 1;
            D1(:,i) = video2patches_fast(cube1,B1,B1,1,1);
        end
    end
 
end