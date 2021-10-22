
%% Simulink Program Hyper-Parameters
T = 0.05;         % Sample Time
NOISE = 0.05;     % Noise

%% MPC Hyper-Paramters


x0 = [1;1;2;3];


% Paramters for Cost function
% Z_r as an inpu

Q = [Qb 0 0 0;
      0 Qa 0 0;
      0 0 Qc 0;
      0 0 0 0];
  
P = Q;
      
R = alpha*eye(1);


    

%% Get State Space
Plant_Model % Get Exact Continuous Plant [A,B,C,D]

A_Aug = C*A*inv(C);
%A_Aug(4,3) = 0;

B_Aug = [0;B(4,1);0;B(3,1)];

C_Aug = eye(4);

D_Aug = [0;0;0;0];

%Choose 1 discretisation
[Aq, Bq, Cq,  Dq] = Plant_Model_ZOH(T,A_Aug,B_Aug,C_Aug,D_Aug); % ZOH exact 
%[Aq, Bq, Cq,  Dq] = Plant_Model_F_Euler(T) % Forward Euler Approximate
%[Aq, Bq, Cq,  Dq] = RLS_Model(T) % From Recursive least sqaures

% %OR Use Model, comment out if you want to use exact model
% % ONly Aq and Bq matter
% Aq = A_discrete_model;
% Bq = B_discrete_model
% Cq = C_discrete_model;
% Dq = D_discrete_model;
% 
% % Augment them to suit correct state space
% Aq = Cq*Aq*inv(Cq);
% Bq = [0;Bq(4,1);0;Bq(3,1)];
% Cq = eye(4);
% Dq = [0;0;0;0];


%% Get Matrices for QuadProg
PHI = cell(N,1); % Each cell within the cell will be 4x4
GAMMA = cell(N,N); % Each cell within the cell will be 1x1


for i=1:N
    PHI{i} = Aq^i; 
    for j=1:N
        if j <= (i-1)
           GAMMA{j,i} = zeros(size(Bq));
        else
            GAMMA{j,i} = (Aq^(j-i))*Bq; 
        end
        
    end
end

PHI = cell2mat(PHI);  % correct
GAMMA = cell2mat(GAMMA);   % correct


OMEGA = [];
PSI = [];

for i=1:N
    PSI = blkdiag(PSI,R);
    if i==N;
        OMEGA = blkdiag(OMEGA,P);
        break
    end
    OMEGA = blkdiag(OMEGA,Q);   
end


%PSI = diag(R,N);


G = 2*(PSI + transpose(GAMMA)*OMEGA*GAMMA);
F = 2*transpose(GAMMA)*OMEGA*PHI;

% G needs to be symettric, It most likely not be due to computational
% roundings

G=(G+G')/2;


% Constraints
Mi = [0 0 0 0; 0 0 0 0; 1 0 0 0;
     -1 0 0 0; 0 1 0 0;0 -1 0 0];
Ei = [1; -1; 0; 0; 0; 0];
fi = [2500;2500;0.2;0.2;0.4;0.4];
Gn = Mi;
h = zeros(size(fi));

D_ = [];
M = [];
E = [];
c = [];

for i=1:N+1 
    if i==1
        fprintf("i==1\n")
        D_ = Mi;
        c = fi;
        %Same as usual
        E = blkdiag(E,Ei);
    elseif i==N+1
        fprintf("i==N\n")
        c = [c;h]
        % Same as usual
        M = blkdiag(M,Gn);
        D_ = [D_;zeros(size(Mi))];
    else
        fprintf("i!=N i!=1\n")
        % Usual
        D_ = [D_;zeros(size(Mi))];
        c = [c;fi];
        M = blkdiag(M,Mi);
        E = blkdiag(E,Ei);
    end
end

M = [zeros(size(Mi,1),size(Mi,2)*N); M];
E = [E; zeros(size(Ei,1),size(Ei,2)*N)];

J = M*GAMMA + E;
W = -M*PHI -D_;


% Assign to Quadprog formulation
H = G; 
f = F*x0;
A = J;
b = c+W*x0;
Aeq = [];
beq= [];
lb = [];
ub = [];

%  quadprog(H,f,A,b,Aeq,beq,lb,ub)

u = quadprog(H,f,A,b)



