function [EoF, M_opt] = ana_two_qubits(Rho, safe_0)
% anlaytic formula of two-qubits state ('Wootters') N=2*2 by concurrence
% Input: density matrix Rho
% Output: Entanglement of formation EoF
%         Optimal cardinality M_opt
% ignore the eigen-values smaller than safe_0

Q = [0 0 0 -1; 0 0 1 0; 0 1 0 0;-1 0 0 0]; % spin-flip transformation
Flip = Q*conj(Rho)*Q;
eig_value = eig(Rho*Flip);
eig_value(abs(eig_value)<safe_0) = 0;
lambda = sort(sqrt(eig_value),'descend'); % square roots of e-values of the non-Hermitian matrix rho*Flip
concurrence = max(0, lambda(1)-lambda(2)-lambda(3)-lambda(4));% 0<=c<=1
x = 1/2 + sqrt(1-concurrence^2)/2; % 1/2<=x<=1
if x == 1 
    EoF = 0;
else
    EoF = -x*log2(x)-(1-x)*log2(1-x); % binary entropy
end

r = rank(Rho);
if r == 3 || lambda(1)<lambda(2)+lambda(3)
    M_opt = 4;
else
    M_opt = r;
end

end
