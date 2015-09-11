function [ state ] = resizeD(state, scale, grow, shrink)
    D = state.D; S = state.S; Z = state.Z;
    Pi = state.Pi; gs = state.gs; geps = state.geps;
    
    [P,K] = size(D);
    B = sqrt(P);
    Pout = (scale*B)^2;
    Dout = [];
    PiOut = [];
    for k=1:K
        dk = reshape(D(:,k), [B,B]);
        if grow && scale < 1
            newD = video2patches_fast(dk, B/2,B/2, 1, 1);
            Dout = [Dout, newD];
            PiOut = [PiOut, Pi(k)*ones(1,size(newD,2))];
        else
            Dout(:,k) = reshape(imresize(dk, scale), Pout,1);
            PiOut(k) = Pi(k);
        end
    end
    PiOut/sum(PiOut);
    if shrink && scale > 1
        q = quantile(state.Pi,20);
        Dout(:,state.Pi <= q(16)) = [];
        state.Pi(state.Pi <= q(16)) = [];
        PiOut = state.Pi;
        disp([K,length(PiOut)]);
    end
    state.D = Dout;
    state.Pi = PiOut;
end

