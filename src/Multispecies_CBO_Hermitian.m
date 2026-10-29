function [EoF_opt,cardinality_opt,Objective] = Multispecies_CBO_Hermitian(rho,N_a,N_b,M_max,method,top_k,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0)

%% multispecies CBO Hermitian method

% dH_j = -lambda*(H_j-barH)dt + sigma*dW_k

% consensus barH is always Hermitian
% method = 'weighted_average'
% (1) barH = sum H_j w_j / sum w_j
%     w_j = exp(-beta*E(H_j))
% method = 'minimizer' (2) barH = argmin E(H_j)

% bipartite quantum state 
% N_a, N_b: dimension of subsytems,  N=N_a*N_b
% rho: N*N density matrix 

%% allow information exchange of different cardinality M
% different M from range [rank(rho), M_max] as a whole
% lift dimension by zero padding
% dimension reduction from M to k
% if top_k = 'eigenvalue': reduce dimension by taking largest k absolute eigenvalues
% if tpo_k = 'norm' : reduce dimension by taking largest k norm column vectors

% output:
% the minimizer EoF_opt and the associated cardinality_opt
% objective: the trajectory of the entanglement at the consensus for all M 


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




%% initialization
rank_state = rank(rho);
M_min = rank_state; 
M_range = M_min:M_max; 
N_size = length(M_range);
particles = cell(N_size,1);
for m = 1:N_size
    
    M = M_range(m);
    H = zeros(M,M,J);
    for j = 1:J
        % generate random Hermitian matrix by GUE: H_{ij} ~ N(0,2/M), forall i,j=1,...,M (?)
        for m1 = 1:M
            H(m1,m1,j) = sqrt(2/M)*randn;
            for m2 = m1+1:M
                H(m1,m2,j) = sqrt(1/M)*randn+1i*sqrt(1/M)*randn;
                H(m2,m1,j) = H(m1,m2,j)';
            end
        end
    end
    particles{m} = H;
    
end


%% preparation
% reduce and lift dimension for every particle in each species
particle_table = cell(N_size);
for n = 1:N_size
    H = particles{n};
    modulated_H = dim_change(H, M_min, M_max, top_k); % cell of {H_m}, H_m is a m*m*J matrix, m = rank,..., 2N
    for m = 1:N_size
        % store particles changed from dimension n+rank-1 to dimension m+rank-1
        particle_table{n,m} = modulated_H{m};
    end
end
% the associated entanglement after dimension change
% note: lift dimension does not change entanglement
Entanglement_table = cell(N_size);
for n = 1:N_size
    for m = 1:N_size
        modulated_H = particle_table{n,m}; % particles from n->m
        Entangle = zeros(1,J);
        for j = 1:J
            Entangle(j) = Entanglement(rho,N_a,N_b,expm(1i*modulated_H(:,:,j)),safe_0);
        end
        Entanglement_table{n,m} = Entangle;
    end
end

% initialized consensus for each species
% take the contribution from other species into account
consensus = cell(N_size,1);
for m = 1:N_size
    % compute consensus for dimension M
    M = M_range(m);
    barH = zeros(M);
    switch method
        case 'weighted_average'
            % (1) barH = sum H_j w_j / sum w_j
            %     w_j = exp(-beta*E(H_j))
            sum_weights = 0;
            for n = 1:N_size
                modulated_H = particle_table{n,m}; % particles from n->m
                Entangle_n = Entanglement_table{n,m}; % associated entanglement
                gibbs_n = exp(-beta*Entangle_n);
                barH = barH + sum(reshape(repelem(gibbs_n,M,M),M,M,J).*modulated_H,3);
                sum_weights = sum_weights + sum(gibbs_n);
            end
            barH = barH/sum_weights;
        case  'minimizer'
            % (2) barH = argmin E(H_j)
            minimizer = Inf;
            for n = 1:N_size
                modulated_H = particle_table{n,m}; % particles from n->m
                Entangle_n = Entanglement_table{n,m}; % associated entanglement
                [min_n,idx] = min(Entangle_n);
                if min_n < minimizer
                    barH = modulated_H(:,:,idx);
                    minimizer = min_n;
                end
            end
    end
    consensus{m} = barH;
end


%% CBO iteration
tic

Objective = zeros(N_size,Kmax); % record the objective value (Entanglement of formation) at the consensus for all dimensions
% min_Entangle = zeros(N_size,Kmax); % record the minimum entanglement of each iteration

for k = 1:Kmax
    
    % update particles
    for m = 1:N_size
        % for each species, drift towards its consensus
        M = M_range(m);
        H = particles{m};
        barH = consensus{m};
        for j = 1:J            
            % generate Hermitian noise
            noise = randn(M,M);
            W = RealToHerm(noise);
            % Euler-Maruyama
            H(:,:,j) = H(:,:,j) - lambda*dt*(H(:,:,j)-barH) + sigma*sqrt(dt)*(a*(H(:,:,j)-barH)+b).*W;            
        end
        particles{m} = H;
    end
    
    
    % update particle and the associated entanglement tables
    % reduce and lift dimension for every particle in each species
    particle_table = cell(N_size);
    for n = 1:N_size
        H = particles{n};
        modulated_H = dim_change(H, M_min, M_max, top_k); % cell of {H_m}, H_m is a m*m*J matrix, m = rank,..., 2N
        for m = 1:N_size
            % store particles changed from dimension n+rank-1 to dimension m+rank-1
            particle_table{n,m} = modulated_H{m};
        end
    end
    % the associated entanglement after dimension change
    % note: lift dimension does not change entanglement
    Entanglement_table = cell(N_size);
    for n = 1:N_size
        for m = 1:N_size
            modulated_H = particle_table{n,m}; % particles from n->m
            Entangle = zeros(1,J);
            for j = 1:J
                Entangle(j) = Entanglement(rho,N_a,N_b,expm(1i*modulated_H(:,:,j)),safe_0);
            end
            Entanglement_table{n,m} = Entangle;
        end
    end
    
    
    % Update consensus
    % take the contribution from other species into account
    consensus = cell(N_size,1);
    for m = 1:N_size
        % compute consensus for dimension M
        M = M_range(m);
        barH = zeros(M);
        switch method
            case 'weighted_average'
                % (1) barH = sum H_j w_j / sum w_j
                %     w_j = exp(-beta*E(H_j))
                sum_weights = 0;
                for n = 1:N_size
                    modulated_H = particle_table{n,m}; % particles from n->m
                    Entangle_n = Entanglement_table{n,m}; % associated entanglement
                    gibbs_n = exp(-beta*Entangle_n);
                    barH = barH + sum(reshape(repelem(gibbs_n,M,M),M,M,J).*modulated_H,3);
                    sum_weights = sum_weights + sum(gibbs_n);
                end
                barH = barH/sum_weights;
            case  'minimizer'
                % (2) barH = argmin E(H_j)
                minimizer = Inf;
                for n = 1:N_size
                    modulated_H = particle_table{n,m}; % particles from n->m
                    Entangle_n = Entanglement_table{n,m}; % associated entanglement
                    [min_n,idx] = min(Entangle_n);
                    if min_n < minimizer
                        barH = modulated_H(:,:,idx);
                        minimizer = min_n;
                    end
                end
        end
        consensus{m} = barH;
    end
    
    
    % record the results
    % Entanglement of the consensus
    Entangle_con = zeros(N_size,1);
    for m = 1:N_size
        U_con = expm(1i*consensus{m});
        Entangle_con(m) = Entanglement(rho,N_a,N_b,U_con,safe_0);
    end
    Objective(:,k) = Entangle_con;
    % Minimum entanglement ( = objective chosing method (2))
    

end


E_M = min(Objective,[],2);
[EoF_opt,index] = min(E_M);
cardinality_opt = M_range(index);

disp(['Multispecies CBO + Hermitian: ', 'rank = ',num2str(rank_state),' cardinality = ',num2str(cardinality_opt),' EoF = ',num2str(EoF_opt),' CPU time = ', num2str(toc)])


end