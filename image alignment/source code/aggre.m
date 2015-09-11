% aggregation the whole image set not using the sift feature
load('C:\Users\haow889\Documents\MATLAB\data\gap160.mat');
[width,height,n]=size(x);
iternum=n-9;
agre_set=zeros(width,height,iternum);
for i=1:iternum
    agre_set(:,:,i)=aggregation(x(:,:,i:i+9),7);
    fname=sprintf('agre%d_%d',i,i+9);
    imwrite(agre_set(:,:,i),strcat('C:\Users\haow889\Documents\MATLAB\',...
        'figure\agre_images\',fname),'png');
end

