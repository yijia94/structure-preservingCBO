function [Entangle,p,Psi] = Entanglement(Rho,N_a,N_b,Unitary,safe_0)
% entanglement of formation
% given a density matrix Rho, compute the entanglement w.r.t a Unitary transformation U
% N_a, N_b : dimension of subsystem 
% associated decomposition: Rho = sum_k p(k)Psi(k)*Psi(k)^H
% ignore the eigen-values smaller than safe_0

% preconditioning : spectral decomposition of density matrix Rho
[V,D] = eig(Rho);  % Rho=V*D*V^H
[d,ind] = sort(diag(D),'descend');
eig_value = d(abs(d)>=safe_0);
r = length(eig_value); % rank(Rho)
V = V(:,ind(1:r)); % first r columns
R = V*diag(sqrt(eig_value)); % Rho=R*R^H

M = length(Unitary);
U = Unitary(:,1:r); % U^H*U=I_r

W = R*U'; % Rho=W*W^H
P = W'*W; % trace(P)=1, P=U*diag(d)*U'
p = diag(P)'; % sum(p)=1


% compute the entanglement of pure states: Psi
E = zeros(1,M);
Psi = zeros(size(W));
for k = 1:M
    
    if p(k) ~= 0
                 
        psi_k = W(:,k)/sqrt(p(k)); % Rho=Psi*diag(p)*Psi^H
        A = reshape(psi_k,N_b,N_a).'; % transpose (no complex conjugate)
        reduced_density_k = A*A'; % trace B subsystem
        eig_value_k = eig(reduced_density_k);
        eig_value_k(abs(eig_value_k)< safe_0) = 0;
        % von-Neumann entropy E(psi_k)
        E(k) = sum(-eig_value_k(eig_value_k~=0).*log2(eig_value_k(eig_value_k~=0))); 
        Psi(:,k) = psi_k;
        
    end
    
end


% the entanglement of mixed states: Rho
Entangle = sum(p.*E);


end








