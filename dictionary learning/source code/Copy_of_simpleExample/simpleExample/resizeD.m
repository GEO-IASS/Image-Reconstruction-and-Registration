function [ Dout ] = resizeD(Din, scale, grow)
    [P,K] = size(Din);
    B = sqrt(P);
    Pout = (scale*B)^2;
    Dout = [];
    for k=1:K
        dk = reshape(Din(:,k), [B,B]);
        if grow && scale < 1
            Dout = [Dout, video2patches_fast(dk, B/2,B/2, 1, 1)];
        else
            Dout(:,k) = reshape(imresize(dk, scale), Pout,1);
        end
    end
end

