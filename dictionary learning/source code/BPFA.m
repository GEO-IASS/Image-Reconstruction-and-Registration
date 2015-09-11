function [state, ss] = BPFA(X,K, Phi, state, ss, opts)
    if nargin<6, opts = struct(); end
    if ~isfield(opts,'ss'), opts.ss = false; end
    if ~isfield(opts,'global'), opts.global = true; end
    if ~isfield(opts,'local'), opts.local = true; end
    if ~isfield(opts, 'init'), opts.init = false; end

    [P,N] = size(X);
    if nargin<4 || isempty(state)
        state.D = randn(P,K);
        state.S = sprandn(N,K,0.2);
        state.Z = full(state.S~=0);
        state.Pi = mean(state.Z);
        state.gs = 1;
        state.geps = 1;
    end
    if nargin<5 || isempty(ss);
        for k=1:K, ss.d_sig{k} = P*(K-k+1); ss.d_mu{k} = 0; end
        ss.pi_a = 1; ss.pi_b = N/4;
        ss.geps_a = 1e-6; ss.geps_b = 1e-6;
        ss.gs_a = 1e-6; ss.gs_b = 1e-6;
    end
    if opts.init, return; end
    
    D = state.D; S = state.S; Z = state.Z;
    Pi = state.Pi; gs = state.gs; geps = state.geps;
    
    R = Phi.*X - Phi.*(D*(Z.*S)');
    resUp = @(updown, R,D,S,Z,Phi,k)... 
        R(:,Z(:,k)) + updown*sparse_mult(Phi(:,Z(:,k)),D(:,k),S(Z(:,k),k));

    if opts.global
        for k=1:K
            R(:,Z(:,k)) = resUp(1, R,D,S,Z,Phi,k);

            sigt = ss.d_sig{k} + geps*Phi(:,Z(:,k))*(S(Z(:,k),k).^2);
            sig_Dk = 1./sigt;        
            mut = ss.d_mu{k} + R(:,Z(:,k))*S(Z(:,k),k);
            mu_D = geps*sig_Dk.*mut;        
            D(:,k) = mu_D + sqrt(sig_Dk).*randn(P,1);

            R(:,Z(:,k)) = resUp(-1, R,D,S,Z,Phi,k);
            if opts.ss, ss.d_sig{k} = sigt; ss.d_mu{k} = mut; end
        end
    end

    if opts.local
        for k=1:K
            R(:,Z(:,k)) = resUp(1, R,D,S,Z,Phi,k);

            DtD = D(:,k)'.^2*Phi;
            DtR = D(:,k)'*R;
            sk = full(S(:,k));
            sk(~Z(:,k)) = randn(N - nnz(Z(:,k)),1)./sqrt(gs);
            probz = -geps*(sk.^2.*DtD'/2 - sk.*DtR');
            probz = exp(probz)*Pi(k);
            Z(:,k) = rand(N,1) > (1-Pi(k))./(probz + 1-Pi(k));

            R(:,Z(:,k)) = resUp(-1, R,D,S,Z,Phi,k);
        end

        for k=1:K
            R(:,Z(:,k)) = resUp(1, R,D,S,Z,Phi,k);

            numz = nnz(Z(:,k));
            DtD = D(:,k)'.^2*Phi(:,Z(:,k));
            DtR = D(:,k)'*R(:,Z(:,k));
            sig_s = 1./(gs + geps*DtD');
            mu_s = geps*sig_s.*DtR';
            S(:,k) = sparse(find(Z(:,k)),1,  mu_s + sqrt(sig_s).*randn(numz,1),   N,1);

            R(:,Z(:,k)) = resUp(-1, R,D,S,Z,Phi,k);
        end
    end
    
    if opts.global
        a = ss.pi_a + sum(Z,1);
        b = ss.pi_b + N - sum(Z,1);
        state.Pi = betarnd(a,b);

        c = ss.geps_a + sum(Phi(:))/2;
        d = ss.geps_b + sum(R(:).^2)/2;
        state.geps = gamrnd(c,1./d);

        e = ss.gs_a + numel(Z)/2;
        f = ss.gs_b + sum(S(:).^2)/2 + (numel(Z)-nnz(Z))/sqrt(gs)/2;
        state.gs = gamrnd(e,1./f);
        %g...,h,i,j,k,l,m,n,o,p... next time won't you sing with me?
        
        if opts.ss
            ss.pi_a = ss.pi_a + a; ss.pi_b = ss.pi_b + b;
            ss.geps_a = ss.geps_a + c; ss.geps_b = ss.geps_b + d;
            ss.gs_a = ss.gs_a + e; ss.gs_b = ss.gs_b + f;
        end
    end
    
    state.D = D; state.S = sparsify(S); state.Z = Z;
end

function DSsparse = sparse_mult(Yflag,D,S)
    [P,N]=size(Yflag);
    DSsparse = sparse(1:P,1:P,double(D(:,1)))*Yflag*sparse(1:N,1:N,double(S(:,1)));
    for k=2:size(D,2)
        DSsparse = DSsparse + sparse(1:P,1:P,double(D(:,k)))*Yflag*sparse(1:N,1:N,double(S(:,k)));
    end
end

function S = sparsify(S)
    [N,K] = size(S);
    [i,j,s] = find(S);
    S = sparse(i,j,s,N,K);
end
