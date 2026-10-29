function [EoF_opt,cardinality_opt,EoF] = M_wrapper(algo,rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0)
%% CBO methods for computing the Engtanglement of Formation with variant cardinality M
% algo: Hermitian preserving, Unitary preserving, Hermitian projection, Unitary projection
% bipartite quantum state 
% N_a, N_b: dimension of subsytems,  N=N_a*N_b
% rho: N*N density matrix 
% theoretical range of cardinality M in [rank(rho), N^2]
% computational range M in [rank(rho), M_max]; pratical cutoff with  M_max = 2*N

% output:
% the minimizer EoF_opt and the associated cardinality_opt
% EoF: the minimizers found for each M in range

rank_state = rank(rho);
M_range = rank_state:1:M_max; 
N_size = length(M_range);
EoF = zeros(N_size,1);

tic
for n = 1:N_size
    
    M = M_range(n); % size of realization
    
    switch algo
        case 'Hermitian'
            Objective = CBO_Hermitian(rho,N_a,N_b,method,M,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
            
        case 'Unitary' % in this case, method is fixed as 'minimizer' 
            Objective = CBO_Unitary(rho,N_a,N_b,M,Kmax,lambda,sigma,a,b,J,dt,safe_0);
            
        case 'Hermitian projection'
            Objective = CBO_Hermitian_projection(rho,N_a,N_b,method,M,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
            
        case 'Unitary projection'
            Objective = CBO_Unitary_projection(rho,N_a,N_b,method,M,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    end
    
    EoF(n) = min(Objective);
end

[EoF_opt, index] = min(EoF);
cardinality_opt = M_range(index);

disp(['CBO + ', algo,': rank = ',num2str(rank_state),' cardinality = ',num2str(cardinality_opt),' EoF = ',num2str(EoF_opt),' CPU time = ', num2str(toc)])


end
