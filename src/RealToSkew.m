function H = RealToSkew(Y)
% Convert a MxM real matrix Y to a SkewHermitian matrix H
% The upper triangular part of Y corresponds to the real part of H
% The lower triangular ... imaginary part ...
M = size(Y,1);
H = zeros(M,M);
for m1 = 1:M
    H(m1,m1) = 1i*Y(m1,m1); % pure imaginary
    for m2 = m1+1:M
        H(m1,m2) = Y(m1,m2)+1i*Y(m2,m1);
        H(m2,m1) = -Y(m1,m2)+1i*Y(m2,m1);
    end
end
end