function result = myresize(vector, num,B,d)
    slicenum = length(vector)/B^2;
    cube1 = patches2video_fast(vector,B,B,slicenum,d,d);
    B2 = B * num;
    cube2 = zeros(B2,B2,slicenum);
    cube2(1:B2/2, 1:B2/2, :)= repmat(cube1(1,1,:),[num num 1]);
    cube2(1:B2/2, B2/2+1:B2, :) = repmat(cube1(1,2,:),[num num 1]);
    cube2(B2/2+1:B2, 1:B2/2, :) = repmat(cube1(2,1,:),[num num 1]);
    cube2(B2/2+1:B2, B2/2+1:B2,:) = repmat(cube1(2,2,:),[num num 1]);
    
    result = video2patches_fast(cube2,B2,B2,d,d);
end