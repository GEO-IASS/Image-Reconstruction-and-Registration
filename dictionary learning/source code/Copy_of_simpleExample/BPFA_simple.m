function [state] = BPFA_simple(X,K, Phi, state)
% Draws a sample from the Beta-process factor analysis model given data X,
% number of dictionary elements K, Binary inpainting mask Phi, and prior state
% (D: dictionary matrix, S: real weight matrix, Z:binary weight-mask matrix, Pi: probability
% vector of dictionary element use, gs: Gamma_S the inverse variance of the
% entries in S, geps: the inverse variance of the noise epsilon).
% 
% Recall, we want to find D,S,Z s.t. X = D(S.*Z) + E
% Iteratively calling this function will eventually provide accurate D,S,Z and the
% inverse variance of E. These are returned in the state variable.

    [P,N] = size(X);
    if isempty(state)
        state.D = randn(P,K);
        state.S = sprandn(N,K,0.2);
        state.Z = full(state.S~=0);
        state.Pi = mean(state.Z);
        state.gs = ones(1,K);
        state.geps = 1;
    end
    % for convenience
    D = state.D; S = state.S; Z = state.Z;
    Pi = state.Pi; gs = state.gs; geps = state.geps;
    
    % hyperparameters
    pi_a = 1;
    pi_b = N/4;
    geps_a = 1e-6;
    geps_b = 1e-6;
    gs_a = 1e-6;
    gs_b = 1e-6;

    % compute residual
    R = Phi.*X - Phi.*(D*(Z.*S)');
    
    % draw each column of D
    for k=1:K
        R(:,Z(:,k)) = R(:,Z(:,k)) + Phi(:,Z(:,k)).*(D(:,k)*S(Z(:,k),k)');

        sig_Dk = 1./(P + geps*Phi(:,Z(:,k))*(S(Z(:,k),k).^2));        
        mu_D = geps*sig_Dk.*( R(:,Z(:,k))*S(Z(:,k),k) );        
        D(:,k) = mu_D + sqrt(sig_Dk).*randn(P,1);

        R(:,Z(:,k)) = R(:,Z(:,k)) - Phi(:,Z(:,k)).*(D(:,k)*S(Z(:,k),k)');
    end

    % draw each row of Z (we've assumed the entries within a row are independent)
    % this is faster since K << N
    for k=1:K
        R(:,Z(:,k)) = R(:,Z(:,k)) + Phi(:,Z(:,k)).*(D(:,k)*S(Z(:,k),k)');

        DtD = D(:,k)'.^2*Phi;
        DtR = D(:,k)'*R;
        sk = full(S(:,k));
        sk(~Z(:,k)) = randn(N - nnz(Z(:,k)),1)./sqrt(gs(k));
        probz = -geps*(sk.^2.*DtD'/2 - sk.*DtR');
        probz = exp(probz)*Pi(k);
        Z(:,k) = rand(N,1) > (1-Pi(k))./(probz + 1-Pi(k));

        R(:,Z(:,k)) = R(:,Z(:,k)) - Phi(:,Z(:,k)).*(D(:,k)*S(Z(:,k),k)');
    end

    % draw each row of S (we've assumed the entries within a row are independent)
    % this is faster since K << N
    for k=1:K
        R(:,Z(:,k)) = R(:,Z(:,k)) + Phi(:,Z(:,k)).*(D(:,k)*S(Z(:,k),k)');

        DtD = D(:,k)'.^2*Phi(:,Z(:,k));
        DtR = D(:,k)'*R(:,Z(:,k));
        sig_s = 1./(gs(k) + geps*DtD');
        mu_s = geps*sig_s.*DtR';
        S(:,k) = sparse(find(Z(:,k)),1,  mu_s + sqrt(sig_s).*randn(nnz(Z(:,k)),1), N,1);

        R(:,Z(:,k)) = R(:,Z(:,k)) - Phi(:,Z(:,k)).*(D(:,k)*S(Z(:,k),k)');
    end

    % draw Pi
    a = pi_a + sum(Z,1);
    b = pi_b + N - sum(Z,1);
    state.Pi = betarnd(a,b);

    % draw gamma_epsilon
    a = geps_a + nnz(Phi)/2;
    b = geps_b + sum(R(:).^2)/2;
    state.geps = gamrnd(a,1./b);

    % draw gamma_S
%     a = gs_a + numel(Z)/2;
%     b = gs_b + sum(S(:).^2)/2 + (numel(Z)-nnz(Z))/sqrt(gs)/2;
%     state.gs = gamrnd(a,1./b);
    a = gs_a + N/2;    
    b = gs_b + sum(S.^2,1)/2 + (N-sum(Z,1))./sqrt(gs)/2;
    state.gs = gamrnd(a,1./b);
    
    state.D = D;
    state.S = sparsify(S);
    state.Z = Z;    
end

function S = sparsify(S)
    [N,K] = size(S);
    [i,j,s] = find(S);
    S = sparse(i,j,s,N,K);
end
