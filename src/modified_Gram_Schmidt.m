function [Q,R] = modified_Gram_Schmidt(A)
% QR decomposition A = QR 
% where the column vectors of Q is orthonormal, R is upper triangle

V = A;
M = length(A);
Q = zeros(size(A));
R = zeros(size(A));
for j = 1:M
    
    R(j,j) = norm(V(:,j));
    q = V(:,j)/R(j,j);
    for k = j+1:M
       R(j,k) = q'*V(:,k);
       V(:,k) = V(:,k)-R(j,k)*q;
    end
    Q(:,j) = q;

end

end