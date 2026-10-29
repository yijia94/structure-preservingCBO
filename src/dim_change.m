function  modulated_particles = dim_change(H, M_min, M_max, top_k)
% change M*M*J matrix H into matrices {H_k}, k = M_min,..., M_max
% where H_k is a k*k*J matrix
% lift dimension by zero padding
% if top_k = 'eigenvalue': reduce dimension by taking largest k absolute eigenvalues
% if tpo_k = 'norm' : reduce dimension by taking largest k norm column vectors

if M_min > M_max
    disp('Warning: M_min should be no greater than M_max!')
end

[M, ~, J] = size(H);
M_size = M_max - M_min + 1;
modulated_particles = cell(M_size,1);


%% reduce dimension M -> k when k < M
for k = M_min:min(M-1,M_max)
    
    H_new = zeros(k,k,J);
    for j = 1:J
        H0 = H(:,:,j);
        switch top_k
            case 'eigenvalue'
                % take largest k absolute eigenvalues
                [V, lambda] = eig(H0, 'vector');
                [~, idx] = sort(abs(lambda), 'descend');
                V_sorted = V(:, idx); % sorted eigen-vectors, V_sorted(:,1:k)'*V_sorted(:,1:k) = I_k
                V_top = V_sorted(1:k, 1:k); % V_top no longer orthogonal
                lambda_top = lambda(idx(1:k));
                H_new(:,:,j) = V_top*diag(lambda_top)*V_top';
            case 'norm' 
                % take largest k norm vectors
                lambda = vecnorm(H0);
                [~, idx] = sort(lambda, 'descend');
                H_sorted = H0(idx,idx); % sorted column vectors
                H_new(:,:,j) = H_sorted(1:k,1:k); % top k submatrix
        end      
    end
    modulated_particles{k-M_min+1} = H_new;
    
end



%% stay H when k = M
if M_min<= M && M<= M_max
    modulated_particles{M-M_min+1} = H;
end

%% lift dimension from M to k when k > M
% fill with 0
for k = (M+1):M_max
    modulated_particles{k-M_min+1} = [H zeros(M,k-M,J); zeros(k-M,k,J)];
end

end