function [Objective, min_Entangle] = CBO_Hermitian(rho,N_a,N_b,method,M,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0)
%% CBO algorithms
% Hermitian preserving CBO
% dH_j = -lambda*(H_j-barH)dt + sigma*dW_k

% consensus barH is always Hermitian
% method = 'weighted_average'
% (1) barH = sum H_j w_j / sum w_j
%     w_j = exp(-beta*E(H_j))

% method = 'minimizer' (2) barH = argmin E(H_j)

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
    Entangle(j) = Entanglement(rho,N_a,N_b,expm(1i*H(:,:,j)),safe_0);
end

% initialized consensus
switch method
    case 'weighted_average'
        % (1) barH = sum H_j w_j / sum w_j
        %     w_j = exp(-beta*E(H_j))
        gibbs = exp(-beta*Entangle);
        consensus = sum(reshape(repelem(gibbs,M,M),M,M,J).*H,3)/sum(gibbs);
    case  'minimizer'
        % (2) barH = argmin E(H_j)
        consensus = H(:,:,Entangle == min(Entangle));
end


%% CBO iteration

Objective = zeros(1,Kmax); % record the objective value (Entanglement of formation) at the consensus
min_Entangle = zeros(1,Kmax); % record the minimum entanglement of each iteration

for k = 1:Kmax
    
    for j = 1:J
        
        Hj = H(:,:,j);
        
        % generate Hermitian noise
        noise = randn(M,M);
        W = RealToHerm(noise);
        
        % Euler-Maruyama
        Hj_new = Hj - lambda*dt*(Hj-consensus) + sigma*sqrt(dt)*(a*(Hj-consensus)+b).*W;
        
        % U_j = exp(iH_j) is unitary
        Uj_new = expm(1i*Hj_new);
        
        % compute the entanglement of density matrix w.r.t. each Unitary transformation
        Entangle(j) = Entanglement(rho,N_a,N_b,Uj_new,safe_0);
        H(:,:,j) = Hj_new;
        
    end
    
    % Update consensus
    switch method
        case 'weighted_average'
            % (1) barH = sum H_j w_j / sum w_j
            %     w_j = exp(-beta*E(H_j))
            gibbs = exp(-beta*Entangle);
            consensus = sum(reshape(repelem(gibbs,M,M),M,M,J).*H,3)/sum(gibbs);
        case  'minimizer'
            % (2) barH = argmin E(H_j)
            consensus = H(:,:,Entangle == min(Entangle));
    end
    
    % record the results
    % Entanglement of the consensus
    U_con = expm(1i*consensus);
    Entangle_con = Entanglement(rho,N_a,N_b,U_con,safe_0);
    Objective(k) = Entangle_con;
    % Minimum entanglement ( = objective chosing method (2))
    min_Entangle(k) = min(Entangle);
    
  
end


end

