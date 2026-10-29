function [EoF_SA,cardinality_opt,optimal_unitary] = SA(Rho,N_a,N_b,M_max,chi0,chi_end,alpha,I_change,I_mat,safe_0)
% Simulated Annealing method for computing the Entanglement of Formation
% for a bipartite quantum state Rho
% N_a, N_b: dimension of subsytems,  N=N_a*N_b
% Rho: N*N density matrix 
% theoretical range of cardinality M in [rank(Rho), N^2]
% computational range M in [rank(Rho), M_max]; pratical cutoff with M_max = 2*N

% SA parameters
% chi0: initial angle
% chi_end: lower bound of angle
% alpha: angle reduction coefficent
% I_change:  angle fixed times
% I_mat: repeating times
% safe_0: tolerance

% output:
% the minimizer EoF_SA 
% the associated unitary matrix optimal_unitary attains EoF_SA and its dimension cardinality_opt


N = N_a*N_b;
rank_state = rank(Rho);
M_range = rank_state:1:M_max; 
N_size = length(M_range);
EoF = zeros(N_size,I_mat);
unitary = zeros(2*N,2*N,I_mat,N_size);

tic
for n = 1:N_size
    
    M = M_range(n); % size of realization
    
    for ii_repeat = 1:I_mat
        
        % initialization of random unitary rotation matrix
        phase = 2*pi*rand(M,1);
        V = diag(exp(1i*phase));
        Entangle = Entanglement(Rho,N_a,N_b,V,safe_0);
        
        chi = chi0;
        while chi >= chi_end
            
            for ii = 1:I_change
                
                % generate random Hermitian matrix by GUE
                H = zeros(M);
                for m1 = 1:M
                    H(m1,m1) = sqrt(2/M)*randn;
                    for m2 = m1+1:M
                        H(m1,m2) = sqrt(1/M)*randn+1i*sqrt(1/M)*randn;
                        H(m2,m1) = H(m1,m2)';
                    end
                end
                % eigenvalue decomposition of H
                [eig_vector,eig_value] = eig(H);
                dV = eig_vector*diag(exp(diag(eig_value)*1i*chi))*eig_vector'; % dV=exp(i chi H) is unitary
                V_new = V*dV;
                Entangle_new = Entanglement(Rho,N_a,N_b,V_new,safe_0);
                
                if Entangle_new < Entangle % accept the move if EoF decreases
                    Entangle = Entangle_new;
                    V = V_new;
                end
                
            end
            
            chi = alpha*chi; % annealing
            
        end
        
        EoF(n,ii_repeat) = Entangle;
        unitary(1:M,1:M,ii_repeat,n) = V;
    end
    
end

[E_M,index_rol] = min(EoF,[],2);
[EoF_SA,index] = min(E_M);
cardinality_opt = M_range(index);
optimal_unitary = unitary(1:cardinality_opt,1:cardinality_opt,index_rol(index),index);

disp(['SA: ',' rank = ',num2str(rank_state),' cardinality = ',num2str(cardinality_opt),' EoF = ',num2str(EoF_SA),' CPU time = ', num2str(toc)])

end








