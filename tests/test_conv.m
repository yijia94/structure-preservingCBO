% Computating Entanglement of Formation for benchmark quantum states 
% test of convergence w.r.t particle number J

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.5,...
    'defaultlinelinewidth',3,'defaultpatchlinewidth',.7)

close all
clear

safe_0 = 1e-10; % tolerance

%%% Examples
% Horodecki state 2*2 (rank 2 state)
N_a = 2; N_b = 2; N = N_a*N_b;
q = 1/4; % 0<q<1  (q=0,1 rank 1)
a = 3/4; % 0<a<1
Rho = [q*a^2 0 0 q*a*sqrt(1-a^2);0 (1-q)*(1-a^2) (1-q)*a*sqrt(1-a^2) 0; 0 (1-q)*a*sqrt(1-a^2) (1-q)*a^2 0; q*a*sqrt(1-a^2) 0 0 q*(1-a^2)];

%% Werner state 2*2 (rank 4 state)
% N_a = 2; N_b = 2; N = N_a*N_b;
% F = 0.7; % 0.5<=F<1 (F=1 rank 1)
% Rho = F*[0 0 0 0; 0 1/2 -1/2 0; 0 -1/2 1/2 0;0 0 0 0] + (1-F)/3*[1 0 0 0; 0 1/2 1/2 0; 0 1/2 1/2 0; 0 0 0 1];

%% random state 2*2
% N_a = 2; N_b = 2; N = N_a*N_b;
% A = randn(N);
% B = A*A'; % Hermitian
% Rho = B/trace(B); % trace 1


% ------------- analytic EoF for 2*2 states by Wootters formula--------------
[EoF_real, M_opt] = ana_two_qubits(Rho, safe_0);
disp(['real entanglement value = ',num2str(EoF_real)])


%% Horodecki state 2*4 (rank 5 entangled PPT state)
% N_a = 2; N_b = 4; N = N_a*N_b;
% b = 0.4; % 0<=b<=1 (b=0 rank=1)
% C = [(1+b)/2 0 0 sqrt(1-b^2)/2; 0 b 0 0;0 0 b 0;sqrt(1-b^2)/2 0 0 (1+b)/2];
% Rho = [b*eye(4) b*diag([1 1 1],1); b*diag([1 1 1],-1) C]/(7*b+1);
% if b == 1 || b == 0
%     EoF_real = 0; % (when b=1, 0 separable)
% end
% no analytical result for 0<b<1


%% Isotropic state N=3*3 (rank 9 state)
% N_a = 3; N_b=3; N = N_a*N_b;
% F = 0.6; % 0<=F<=1 ( F = 0 rank = 8 hard, F=1, rank = 1 else rank = 9)
% psi = 1/sqrt(3)*[1 0 0 0 1 0 0 0 1];
% temp = psi'*psi;
% Rho = (1-F)/8*(eye(9)-temp)+F*temp;
% 
% % %------------- analytic by Terhal and Vollbrecht¡¯s formula --------------
% if F > 8/9
%     EoF_real = 3*(F-1)+log2(3);
% elseif F > 1/3
%     gamma = (2-F+sqrt(2*F*(1-F))*2)/3;
%     H = -gamma*log2(gamma)-(1-gamma)*log2(1-gamma);
%     EoF_real = H+(1-gamma);
% else
%     EoF_real = 0;
% end
% disp(['real entanglement value = ',num2str(EoF_real)])








%% test convergence
J = [12,25,50,100,200,400];
rank_state = rank(Rho); % rank
M_max = 2*N; % r<=M<=2N

% CBO parameters
Kmax = 500; % iteration
beta = 200; % inverse temperature >>1
lambda = 1; % drift rate 2lambda> sigma^2
sigma = 0.06; % noise level
a = 1; % cbo noise level
b = 1; % addictive noise level
dt = 0.2; % time step
method = 'weighted_average';



% CBO H
N_J = length(J);
EoF = zeros(N_J,1);
for j = 1:N_J
    EoF(j) = M_wrapper('Hermitian',Rho,N_a,N_b,M_max,method,Kmax,beta,lambda,sigma,a,b,J(j),dt,safe_0);
end
error_H = EoF - EoF_real;
figure(1)
plot(J,error_H)
xlabel('J')
ylabel('Error')


% % multispecies CBO H
top_k = 'norm'; % other choice 'eigenvalue'
EoF2 = zeros(N_J,1);
for j = 1:N_J
    EoF2(j) = Multispecies_CBO_Hermitian(Rho,N_a,N_b,M_max,method,top_k,Kmax,beta,lambda,sigma,a,b,J(j),dt,safe_0);
end
error_multi_H = EoF2 - EoF_real;
figure(2)
plot(J,error_multi_H)
xlabel('J')
ylabel('Error')


