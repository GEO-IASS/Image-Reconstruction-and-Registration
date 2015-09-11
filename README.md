<<<<<<< HEAD
=======
# Image-Reconstruction-and-Registration
>>>>>>> fa092e01f2979bd0c2ddea84fa025dad625fb61a
# First Part
This part is for the dictionary learning. The basic idea is based on the method in reference. Here we tried to speed the original algorithm in a hierarchical way.
#### Goal:
To reconstruct the image with lower pixel fraction and get acceptable psnr value

#### Method:
1. Start from the quarter size of the original image to get the first level dictionary
2. Get the second level dictionary based on former level for the half size image
3. Get the last level dictionary based on former level for the whole size image

In particular, there are different ways to enlarge the dictionary, resample the element coefficients as well as mix enlarge and split procedures.

#### Problem:
In my view, quarter size image has lost a lot of information. Dictionary starting from this level will not work well. The hierarchical idea does speed up the algorithm a lot but does not give us reconstruction of better quality. 

#### Furture Work:
The distributed realization of the original algorithm


# Second Part
This part is the image alignment. The basic idea is based on the method in reference but the result is not as good as expected.

#### Method:
1. Extract and store the sift feature information of each image
2. According to the feature distance figure out the matched image pairs
3. Align the image based on the translation of fixed points and moving points
4. Fourier Burst Accumulation

In particular, Fourier Burst Accumulation is weighted Fourier reconstruction.  I first used FBA directly on the original image set. Then I added registration based on sift features before FBA.

#### Problem:
The result image may have wide black edges or blurred effect. I tried to use different distance definition and distance ratio to further restrict the matched image pairs. There is no quantitative evaluation for good or bad results. Further evaluation is pending. The result images are uploaded on share folder.

#### Future Work:
<<<<<<< HEAD
Try to use different features to figure out the translation of points.Adjust the variables like the distance ratio or distance to get better results.

=======
Try to use different features to figure out the translation of points.
Adjust the variables like the distance ratio or distance to get better results.
>>>>>>> fa092e01f2979bd0c2ddea84fa025dad625fb61a
