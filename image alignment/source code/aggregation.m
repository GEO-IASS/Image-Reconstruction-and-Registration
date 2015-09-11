function aggreim=aggregation(image_set,p_value)
    %takes a series of images as input and outputs the aggregation image
    %p should range from 7 to 20 according to the paper
    
    [width,height,n]=size(image_set);
    weight=zeros(width,height);
    u_p_hat=zeros(width,height,1);
    for i=1:n
        v_i=image_set(:,:,i);
        fv_i=fft2(v_i);
        w_i=abs(fv_i);
        sigma=min(width,height)/50;
        H=fspecial('gaussian',[3 3],sigma);
        smw_i=(imfilter(w_i,H)).^p_value;
        u_p_hat=u_p_hat+smw_i.*fv_i;
        weight=weight+smw_i;
    end
    aggreim=ifft2(u_p_hat./weight);
    
end