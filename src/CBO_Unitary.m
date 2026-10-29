function [min_Entangle] = CBO_Unitary(rho,N_a,N_b,M,Kmax,lambda,sigma,a,b,J,dt,safe_0)
%% CBO algorithms
% Unitary preserving CBO
% dU_j = lambda*(barU U_j^* - U_j barU^*)U_j dt + sigma dW_j U_j

% consensus
% (1) barU = sum U_j w_j / sum w_j (non-Unitary)
%     w_j = exp(-beta*E(U_j))

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
[~,idx] = min(Entangle);
U_con = U(:,:,idx);

%% CBO iteration


% Objective = zeros(1,Kmax); % record the objective value (Entanglement of formation) at the consensus
min_Entangle = zeros(1,Kmax); % record the minimum entanglement of each iteration

for k = 1:Kmax
    
    for j = 1:J
        
        Uj = U(:,:,j);
        
        % generate Skew-Hermitian noise
        noise = randn(M,M);
        W = RealToSkew(noise);
        diffusion = sqrt(dt) * sigma * (a*norm(Uj-U_con)+b) * W; % Skew-Hermitian
        
        % Drift toward consensus
        Uskew = U_con*Uj' - Uj*U_con';
        drift = dt * lambda * Uskew; % Skew-Hermitian
        
        % Euler step
        Uj_new = expm(drift+diffusion)*Uj;
        
        % compute the entanglement of density matrix w.r.t. each Unitary transformation
        Entangle(j) = Entanglement(rho,N_a,N_b,Uj_new,safe_0);
        
        U(:,:,j) = Uj_new;
        
    end
    
    % Update of consensus: method = 'minimizer'
    % barU = argmin E(U_j)
    % Minimum entanglement ( = objective chosing method (2))
    [Entangle_min, idx] = min(Entangle);
    U_con = U(:,:,idx);
    min_Entangle(k) = Entangle_min;

end


end

