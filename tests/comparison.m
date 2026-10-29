%% performance of different optimization methods for EoF

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.5,...
    'defaultlinelinewidth',3,'defaultpatchlinewidth',.7)

close all
clear
tic

parpool(4)


% CBO parameters
beta = 200; % inverse temperature >>1
lambda = 1; % drift rate 2 lambda> sigma^2
sigma = 0.06; % noise level for Hermitian formulation
sigma2 = 0.01; % noise level for Unitary formulation
a = 1; % cbo noise level
b = 1; % addictive noise level
J = 500; % number of particles 
dt = 0.2; % time step
Kmax = 1000; % maximum iteration step
safe_0 = 1e-10; % tolerance
% consensus 
method = 'weighted_average'; % other choice 'minimizer'
% dimension reduction technique for multi-species cbo
top_k = 'norm'; % other choice 'eigenvalue'

% SA parameters
chi0 = 0.3; % initial angle
chi_end = 0.0001; % lower bound of angle
alpha = 2/3; % angle reduction coefficent
I_change = 1000; % angle fixed times
I_mat = 20; % repeating times



%%%  benchmark tests
%% ---------------------ex1-------------------------------------------
N_a = 2; N_b = 2; N = N_a*N_b;
dq = 0.02;
q = 0:dq:1; % 0<q<1
para_a = 3/4; % 0<a<1
M_max = 2*N; % maximum cardinality

N_test = length(q);
real_entangle = zeros(1,N_test); real_cardinality = zeros(1,N_test); rank_state = zeros(1,N_test);
Entangle_SA = zeros(1,N_test); cardinality_SA = zeros(1,N_test);
Entangle_multi_CBO_H = zeros(1,N_test); cardinality_multi_CBO_H = zeros(1,N_test);
Entangle_CBO_H = zeros(1,N_test); cardinality_CBO_H = zeros(1,N_test);
Entangle_CBO_H_proj = zeros(1,N_test); cardinality_CBO_H_proj = zeros(1,N_test);
Entangle_CBO_U = zeros(1,N_test); cardinality_CBO_U = zeros(1,N_test);
Entangle_CBO_U_proj = zeros(1,N_test); cardinality_CBO_U_proj = zeros(1,N_test);

parfor ii_test = 1:N_test
    
    % Horodecki state N=2*2
    Rho = [q(ii_test)*para_a^2 0 0 q(ii_test)*para_a*sqrt(1-para_a^2);0 (1-q(ii_test))*(1-para_a^2) (1-q(ii_test))*para_a*sqrt(1-para_a^2) 0; 0 (1-q(ii_test))*para_a*sqrt(1-para_a^2) (1-q(ii_test))*para_a^2 0; q(ii_test)*para_a*sqrt(1-para_a^2) 0 0 q(ii_test)*(1-para_a^2)];
    rk = rank(Rho);
    rank_state(ii_test) = rk;
    % ------------- analytic EoF for 2*2 states by Wootters formula--------------
    [real_entangle(ii_test), real_cardinality(ii_test)] = ana_two_qubits(Rho, safe_0);
    
    
    [EoF_SA, M_SA] = SA(Rho,N_a,N_b,M_max,chi0,chi_end,alpha,I_change,I_mat,safe_0);
    Entangle_SA(ii_test) = EoF_SA; cardinality_SA(ii_test) = M_SA;
    
    [EoF_multi_CBO_H, M_multi_CBO_H] = Multispecies_CBO_Hermitian(Rho,N_a,N_b,M_max,method,top_k,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_multi_CBO_H(ii_test) = EoF_multi_CBO_H; cardinality_multi_CBO_H(ii_test) = M_multi_CBO_H;
    
    [EoF_CBO_H, M_CBO_H] = M_wrapper('Hermitian',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H(ii_test) = EoF_CBO_H; cardinality_CBO_H(ii_test) = M_CBO_H;
    
    [EoF_CBO_H_proj, M_CBO_H_proj] = M_wrapper('Hermitian projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H_proj(ii_test) = EoF_CBO_H_proj; cardinality_CBO_H_proj(ii_test) = M_CBO_H_proj;
    
    [EoF_CBO_U, M_CBO_U] = M_wrapper('Unitary',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U(ii_test) = EoF_CBO_U; cardinality_CBO_U(ii_test) = M_CBO_U;
    
    [EoF_CBO_U_proj, M_CBO_U_proj] = M_wrapper('Unitary projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U_proj(ii_test) = EoF_CBO_U_proj; cardinality_CBO_U_proj(ii_test) = M_CBO_U_proj;
    
    
end
toc
error_SA = abs(Entangle_SA-real_entangle);
error_multi_CBO_H = abs(Entangle_multi_CBO_H-real_entangle);
error_CBO_H = abs(Entangle_CBO_H-real_entangle);
error_CBO_U = abs(Entangle_CBO_U-real_entangle);
error_CBO_H_proj = abs(Entangle_CBO_H_proj-real_entangle);
error_CBO_U_proj = abs(Entangle_CBO_U_proj-real_entangle);
save result_Horodecki

figure(1)
subplot(1,2,1)
plot(q,Entangle_SA,'v-',q,Entangle_CBO_U,'o--',q,Entangle_CBO_H,'o--',q,Entangle_CBO_U_proj,'o--',q,Entangle_CBO_H_proj,'o--',q,Entangle_multi_CBO_H,'o--',q,real_entangle,'k-')
xlabel('q')
ylabel('EoF')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres','analytic')
title(['Horodecki state 2*2, a =',num2str(para_a)])
subplot(1,2,2)
plot(q,error_SA,'v-',q,error_CBO_U,'o--',q,error_CBO_H,'o--',q,error_CBO_U_proj,'o--',q,error_CBO_H_proj,'o--',q,error_multi_CBO_H,'o--')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres')
xlabel('q')
ylabel('Error')
title(['Horodecki state 2*2, a =',num2str(para_a)])





%% ---------------------ex2-------------------------------------------
N_a = 2; N_b = 2; N = N_a*N_b;
dF = 0.01;
F = 0.5:dF:1; % 0.5<=F<=1
Bell1 = [0 0 0 0; 0 1/2 -1/2 0; 0 -1/2 1/2 0;0 0 0 0];
Bells = [1 0 0 0; 0 1/2 1/2 0; 0 1/2 1/2 0; 0 0 0 1];
M_max = 2*N; % maximum cardinality

N_test = length(F);
real_entangle = zeros(1,N_test); real_cardinality = zeros(1,N_test); rank_state = zeros(1,N_test);
Entangle_SA = zeros(1,N_test); cardinality_SA = zeros(1,N_test);
Entangle_multi_CBO_H = zeros(1,N_test); cardinality_multi_CBO_H = zeros(1,N_test);
Entangle_CBO_H = zeros(1,N_test); cardinality_CBO_H = zeros(1,N_test);
Entangle_CBO_H_proj = zeros(1,N_test); cardinality_CBO_H_proj = zeros(1,N_test);
Entangle_CBO_U = zeros(1,N_test); cardinality_CBO_U = zeros(1,N_test);
Entangle_CBO_U_proj = zeros(1,N_test); cardinality_CBO_U_proj = zeros(1,N_test);

parfor ii_test = 1:N_test
    
    % Werner state 2*2
    Rho = F(ii_test)*Bell1 + (1-F(ii_test))/3*Bells;
    rk = rank(Rho);
    rank_state(ii_test) = rk;
    % ------------- analytic EoF for 2*2 states by Wootters formula--------------
    [real_entangle(ii_test), real_cardinality(ii_test)] = ana_two_qubits(Rho, safe_0);
    
    
    [EoF_SA, M_SA] = SA(Rho,N_a,N_b,M_max,chi0,chi_end,alpha,I_change,I_mat,safe_0);
    Entangle_SA(ii_test) = EoF_SA; cardinality_SA(ii_test) = M_SA;
    
    [EoF_multi_CBO_H, M_multi_CBO_H] = Multispecies_CBO_Hermitian(Rho,N_a,N_b,M_max,method,top_k,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_multi_CBO_H(ii_test) = EoF_multi_CBO_H; cardinality_multi_CBO_H(ii_test) = M_multi_CBO_H;
    
    [EoF_CBO_H, M_CBO_H] = M_wrapper('Hermitian',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H(ii_test) = EoF_CBO_H; cardinality_CBO_H(ii_test) = M_CBO_H;
    
    [EoF_CBO_H_proj, M_CBO_H_proj] = M_wrapper('Hermitian projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H_proj(ii_test) = EoF_CBO_H_proj; cardinality_CBO_H_proj(ii_test) = M_CBO_H_proj;
    
    [EoF_CBO_U, M_CBO_U] = M_wrapper('Unitary',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U(ii_test) = EoF_CBO_U; cardinality_CBO_U(ii_test) = M_CBO_U;
    
    [EoF_CBO_U_proj, M_CBO_U_proj] = M_wrapper('Unitary projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U_proj(ii_test) = EoF_CBO_U_proj; cardinality_CBO_U_proj(ii_test) = M_CBO_U_proj;
    
    
end
toc
error_SA = abs(Entangle_SA-real_entangle);
error_multi_CBO_H = abs(Entangle_multi_CBO_H-real_entangle);
error_CBO_H = abs(Entangle_CBO_H-real_entangle);
error_CBO_U = abs(Entangle_CBO_U-real_entangle);
error_CBO_H_proj = abs(Entangle_CBO_H_proj-real_entangle);
error_CBO_U_proj = abs(Entangle_CBO_U_proj-real_entangle);
save result_Werner

figure(2)
subplot(1,2,1)
plot(F,Entangle_SA,'v-',F,Entangle_CBO_U,'o--',F,Entangle_CBO_H,'o--',F,Entangle_CBO_U_proj,'o--',F,Entangle_CBO_H_proj,'o--',F,Entangle_multi_CBO_H,'o--',F,real_entangle,'k-')
xlabel('F')
ylabel('EoF')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres','analytic')
title('Werner state')
subplot(1,2,2)
plot(F,error_SA,'v-',F,error_CBO_U,'o--',F,error_CBO_H,'o--',F,error_CBO_U_proj,'o--',F,error_CBO_H_proj,'o--',F,error_multi_CBO_H,'o--')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres')
xlabel('F')
ylabel('Error')
title('Werner state')






%% ---------------------ex3-------------------------------------------
N_a = 2; N_b = 4; N = N_a*N_b;
db = 0.02;
para_b = 0:db:1; % 0<b<1
M_max = 2*N; % maximum cardinality

N_test = length(para_b);
rank_state = zeros(1,N_test);
Entangle_SA = zeros(1,N_test); cardinality_SA = zeros(1,N_test);
Entangle_multi_CBO_H = zeros(1,N_test); cardinality_multi_CBO_H = zeros(1,N_test);
Entangle_CBO_H = zeros(1,N_test); cardinality_CBO_H = zeros(1,N_test);
Entangle_CBO_H_proj = zeros(1,N_test); cardinality_CBO_H_proj = zeros(1,N_test);
Entangle_CBO_U = zeros(1,N_test); cardinality_CBO_U = zeros(1,N_test);
Entangle_CBO_U_proj = zeros(1,N_test); cardinality_CBO_U_proj = zeros(1,N_test);

parfor ii_test = 1:N_test
    
    % Horodecki state N=2*4
    bb = para_b(ii_test);
    C = [(1+bb)/2 0 0 sqrt(1-bb^2)/2; 0 bb 0 0;0 0 bb 0; sqrt(1-bb^2)/2 0 0 (1+bb)/2];
    Rho = [bb*eye(4) bb*diag([1 1 1],1); bb*diag([1 1 1],-1) C]/(7*bb+1);
    rk = rank(Rho);
    rank_state(ii_test) = rk;
    
    
    [EoF_SA, M_SA] = SA(Rho,N_a,N_b,M_max,chi0,chi_end,alpha,I_change,I_mat,safe_0);
    Entangle_SA(ii_test) = EoF_SA; cardinality_SA(ii_test) = M_SA;
    
    [EoF_multi_CBO_H, M_multi_CBO_H] = Multispecies_CBO_Hermitian(Rho,N_a,N_b,M_max,method,top_k,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_multi_CBO_H(ii_test) = EoF_multi_CBO_H; cardinality_multi_CBO_H(ii_test) = M_multi_CBO_H;
    
    [EoF_CBO_H, M_CBO_H] = M_wrapper('Hermitian',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H(ii_test) = EoF_CBO_H; cardinality_CBO_H(ii_test) = M_CBO_H;
    
    [EoF_CBO_H_proj, M_CBO_H_proj] = M_wrapper('Hermitian projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H_proj(ii_test) = EoF_CBO_H_proj; cardinality_CBO_H_proj(ii_test) = M_CBO_H_proj;
    
    [EoF_CBO_U, M_CBO_U] = M_wrapper('Unitary',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda2,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U(ii_test) = EoF_CBO_U; cardinality_CBO_U(ii_test) = M_CBO_U;
    
    [EoF_CBO_U_proj, M_CBO_U_proj] = M_wrapper('Unitary projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda2,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U_proj(ii_test) = EoF_CBO_U_proj; cardinality_CBO_U_proj(ii_test) = M_CBO_U_proj;
    
    
end
toc
save result_Horodecki2times4

figure(3)
plot(para_b,Entangle_SA,'v-',para_b,Entangle_CBO_U,'o--',para_b,Entangle_CBO_H,'o--',para_b,Entangle_CBO_U_proj,'o--',para_b,Entangle_CBO_H_proj,'o--',para_b,Entangle_multi_CBO_H,'o--')
xlabel('b')
ylabel('EoF')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres')
title('Horodecki state 2*4')





%% ---------------------ex4-------------------------------------------
N_a = 3; N_b=3; N = N_a*N_b;
psi = 1/sqrt(3)*[1 0 0 0 1 0 0 0 1];
temp = psi'*psi;
dF = 0.02;
F = 0:dF:1; % 0<=F<=1
M_max = 2*N; % maximum cardinality

N_test = length(F);
real_entangle = zeros(1,N_test); rank_state = zeros(1,N_test);
Entangle_SA = zeros(1,N_test); cardinality_SA = zeros(1,N_test);
Entangle_multi_CBO_H = zeros(1,N_test); cardinality_multi_CBO_H = zeros(1,N_test);
Entangle_CBO_H = zeros(1,N_test); cardinality_CBO_H = zeros(1,N_test);
Entangle_CBO_H_proj = zeros(1,N_test); cardinality_CBO_H_proj = zeros(1,N_test);
Entangle_CBO_U = zeros(1,N_test); cardinality_CBO_U = zeros(1,N_test);
Entangle_CBO_U_proj = zeros(1,N_test); cardinality_CBO_U_proj = zeros(1,N_test);

parfor ii_test = 1:N_test
    
    % Isotropic state N=3*3
    Rho = (1-F(ii_test))/8*(eye(9)-temp) + F(ii_test)*temp;
    rk = rank(Rho);
    rank_state(ii_test) = rk;
    
    % analytic by Terhal and Vollbrechts formula
    if F(ii_test) > 8/9
        real_entangle(ii_test) = 3*(F(ii_test)-1)+log2(3);
    elseif F(ii_test) > 1/3
        gamma = (2-F(ii_test)+sqrt(2*F(ii_test)*(1-F(ii_test)))*2)/3;
        H = -gamma*log2(gamma)-(1-gamma)*log2(1-gamma);
        real_entangle(ii_test) = H+(1-gamma);
    else
        real_entangle(ii_test) = 0;
    end
    
    [EoF_SA, M_SA] = SA(Rho,N_a,N_b,M_max,chi0,chi_end,alpha,I_change,I_mat,safe_0);
    Entangle_SA(ii_test) = EoF_SA; cardinality_SA(ii_test) = M_SA;
    
    [EoF_multi_CBO_H, M_multi_CBO_H] = Multispecies_CBO_Hermitian(Rho,N_a,N_b,M_max,method,top_k,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_multi_CBO_H(ii_test) = EoF_multi_CBO_H; cardinality_multi_CBO_H(ii_test) = M_multi_CBO_H;
    
    [EoF_CBO_H, M_CBO_H] = M_wrapper('Hermitian',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H(ii_test) = EoF_CBO_H; cardinality_CBO_H(ii_test) = M_CBO_H;
    
    [EoF_CBO_H_proj, M_CBO_H_proj] = M_wrapper('Hermitian projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J,dt,safe_0);
    Entangle_CBO_H_proj(ii_test) = EoF_CBO_H_proj; cardinality_CBO_H_proj(ii_test) = M_CBO_H_proj;
    
    [EoF_CBO_U, M_CBO_U] = M_wrapper('Unitary',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U(ii_test) = EoF_CBO_U; cardinality_CBO_U(ii_test) = M_CBO_U;
    
    [EoF_CBO_U_proj, M_CBO_U_proj] = M_wrapper('Unitary projection',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma2,a,b,J,dt,safe_0);
    Entangle_CBO_U_proj(ii_test) = EoF_CBO_U_proj; cardinality_CBO_U_proj(ii_test) = M_CBO_U_proj;
    
    
end
toc
error_SA = abs(Entangle_SA-real_entangle);
error_multi_CBO_H = abs(Entangle_multi_CBO_H-real_entangle);
error_CBO_H = abs(Entangle_CBO_H-real_entangle);
error_CBO_U = abs(Entangle_CBO_U-real_entangle);
error_CBO_H_proj = abs(Entangle_CBO_H_proj-real_entangle);
error_CBO_U_proj = abs(Entangle_CBO_U_proj-real_entangle);
save result_Isotropic

figure(4)
subplot(1,2,1)
plot(F,Entangle_SA,'v-',F,Entangle_CBO_U,'o--',F,Entangle_CBO_H,'o--',F,Entangle_CBO_U_proj,'o--',F,Entangle_CBO_H_proj,'o--',F,Entangle_multi_CBO_H,'o--',F,real_entangle,'k-')
xlabel('F')
ylabel('EoF')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres','analytic')
title('Isotropic state')
subplot(1,2,2)
plot(F,error_SA,'v-',F,error_CBO_U,'o--',F,error_CBO_H,'o--',F,error_CBO_U_proj,'o--',F,error_CBO_H_proj,'o--',F,error_multi_CBO_H,'o--')
legend('SA','CBO U pres','CBO H pres','CBO U proj','CBO H proj','Multi-species CBO H pres')
xlabel('F')
ylabel('Error')
title('Isotropic state')



delete(gcp('nocreate'))



