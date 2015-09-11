

function  recoverImage = alignment(siftinfo1,siftinfo2,DisRatio)
% Find SIFT keypoints for Image1 and Image2, the distance threshold would be 
% the disRatio (the angle of feature vectors)
im1 = cell2mat(siftinfo1(1,1));
des1 = cell2mat(siftinfo1(2,1));
loc1 = cell2mat(siftinfo1(3,1));

im2 = cell2mat(siftinfo2(1,1));
des2 = cell2mat(siftinfo2(2,1));
loc2 = cell2mat(siftinfo2(3,1));


% For efficiency in Matlab, it is cheaper to compute dot products between
%  unit vectors rather than Euclidean distances.  Note that the ratio of 
%  angles (acos of dot products of unit vectors) is a close approximation
%  to the ratio of Euclidean distances for small angles.
%
% distRatio: Only keep matches in which the ratio of vector angles from the
%   nearest to second nearest neighbor is less than distRatio.
distRatio = DisRatio;   

% For each descriptor in the first image, select its match to second image.
des2t = des2'; % Precompute matrix transpose
matches=zeros(1,size(des1,1));
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;        % Computes vector of dot products
   [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

   % Check if nearest neighbor has angle less than distRatio times 2nd.
   if (vals(1) < distRatio * vals(2))
      matches(i) = indx(1);
   else
      matches(i) = 0;
   end
end


num = sum(matches > 0);
fixedPointsSet = zeros(num,2);
movingPointsSet = zeros(num,2);
distances = zeros(num,1);
dummy_n = 1;
for i = 1: size(des1,1)
  if (matches(i) > 0)
    point1 = [loc1(i,1) loc1(i,2)];
    point2 = [loc2(matches(i),1) loc2(matches(i),2)];    
    dist = sqrt(sum((point2 - point1).^2));
%     dist = max(abs(point1 - point2)); % use the maximum 
%     dist = sum(abs(point1 - point2)); % use the mahattan distance
    distances(dummy_n,:) = dist;
    fixedPointsSet(dummy_n,:) = [loc1(i,1) loc1(i,2)];
    movingPointsSet(dummy_n,:) = [loc2(matches(i),1) loc2(matches(i),2)];
    dummy_n = dummy_n + 1;
%     if (dist <=10)
%         distances(dummy_n,:) = dist;
%         fixedPointsSet(dummy_n,:) = [loc1(i,1) loc1(i,2)];
%         movingPointsSet(dummy_n,:) = [loc2(matches(i),1) loc2(matches(i),2)];
%         dummy_n = dummy_n + 1;
%     end
  end
end


% recover image
if min(distances) > 10
    recoverImage = im2;
else
    [minimum, idx] = min(distances);
    fixedPoints = fixedPointsSet(idx,:);
    movingPoints = movingPointsSet(idx,:);

    theta = 0;
    tx = fixedPoints(1,1) - movingPoints(1,1);
    ty = fixedPoints(1,2) - movingPoints(1,2);

    tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; tx ty 1]);

% if (size(nonzeros(fixedPoints),1) < 4)
%     recoverImage = im2;
% else
%     tform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
%     tformInv = invert(tform);
%     Tinv = tformInv.T;
%     ss = Tinv(2,1);
%     sc = Tinv(1,1);
%     tx = Tinv(3,1);
%     ty = Tinv(3,2);
%     scale_recovered = sqrt(ss * ss + sc * sc);
%     theta_recovered = atan2(ss,sc) * 180/pi;
%     translation_recovered = [tx ty];
    Roriginal = imref2d(size(im1));
    recoverImage = imwarp(im2, tform, 'OutputView', Roriginal);
% end
end


% to display the final result of the recoverred image
% figure;
% norm(double(im1-im2),'fro')
% norm(double(im1-recovered),'fro')
% imshowpair(im1,recovered)

end

