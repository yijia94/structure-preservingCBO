function [Objective, min_Entangle] = CBO_Unitary_projection(rho,N_a,N_b,method,M,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0)
% % CBO algorithms
% Unitary projection CBO
% dU_j = -lambda*(U_j-barU)dt + sigma*dW
%  U_j non-Unitary,  gram - schmidt to ensure orthogonality

% consensus
% method = 'weighted_average'
% (1) barU = sum U_j w_j / sum w_j (non-Unitary)
%     w_j = exp(-beta*E(U_j))
%     gram - schimdt to ensure orthogonality

% method = 'minimizer' 
% (2) barU = argmin E(U_j)  (Unitary)


% CBO parameters
% beta: inverse temperature >>1
% lambda: drift rate 
% sigma : noise level, must satisfy 2lambda> sigma^2
% a = 1: cbo noise level
% b: addictive noise level
% J:  number of particles
% dt: time step
% Kmax: maximum iteration step
% safe_0: tolerance

% bipartite quantum state 
% N_a, N_b: dimension of subsytems,  N=N_a*N_b
% rho: N*N density matrix 
% M: cardinality (number of pure state decomposition used to create rho)

%% initialization
H = zeros(M,M,J);
U = zeros(M,M,J);
Entangle = zeros(1,J);

for j = 1:J
    % generate random Hermitian matrix by GUE: H_{ij} ~ N(0,2/M), forall i,j=1,...,M 
    for m1 = 1:M
        H(m1,m1,j) = sqrt(2/M)*randn;
        for m2 = m1+1:M
            H(m1,m2,j) = sqrt(1/M)*randn+1i*sqrt(1/M)*randn;
            H(m2,m1,j) = H(m1,m2,j)';
        end
    end
    U(:,:,j) = expm(1i*H(:,:,j));
    Entangle(j) = Entanglement(rho,N_a,N_b,U(:,:,j),safe_0);
end

% initialized consensus
switch method
    case 'weighted_average'
        gibbs = exp(-beta*Entangle);
        consensus = sum(reshape(repelem(gibbs,M,M),M,M,J).*U,3)/sum(gibbs);
        consensus = modified_Gram_Schmidt(consensus);
    case  'minimizer'
        consensus = U(:,:,Entangle == min(Entangle));
end

%% CBO iteration


Objective = zeros(1,Kmax); % record the objective value (Entanglement of formation) at the consensus
min_Entangle = zeros(1,Kmax); % record the minimum entanglement of each iteration

for k = 1:Kmax
    
    
    % Update of particles
    for j = 1:J
        
        Uj = U(:,:,j);
        
        % Euler-Maruyama
        noise = randn(M,M);
        Uj_new = Uj - lambda*dt*(Uj-consensus) + sigma*sqrt(dt)*(a*(Uj-consensus)+b).*noise;
        
        % check linearly independence of the column vectors of Uj_new
        if rank(Uj_new) < M
            disp(['Warning: the ', num2str(j),' -th particle in step ',num2str(k),' is linearly dependent! '])
            Uj_new = U(:,:,j);
        else
            % Schmidt Orthogonal process
            Uj_new = modified_Gram_Schmidt(Uj_new); % Unitary
            Entangle(j) = Entanglement(rho,N_a,N_b,Uj_new,safe_0); % entanglement w.r.t. new Unitary
        end
        U(:,:,j) = Uj_new;
        
    end
    
    % Update consensus
    switch method
        case 'weighted_average'
            % (1) barU = sum U_j w_j / sum w_j
            %   w_j = exp(-beta*E(U_j))
            gibbs = exp(-beta*Entangle);
            consensus = sum(reshape(repelem(gibbs,M,M),M,M,J).*U,3)/sum(gibbs);
            consensus = modified_Gram_Schmidt(consensus);
        case  'minimizer'
            % (2) barU = argmin E(U_j)
            consensus = U(:,:, Entangle==min(Entangle));
    end
    
    
    % record the results
    % Entanglement of the consensus
    Entangle_con = Entanglement(rho,N_a,N_b,consensus,safe_0);
    Objective(k) = Entangle_con;
    % Minimum entanglement ( = objective chosing method (2))
    min_Entangle(k) = min(Entangle);
    

    
end


end