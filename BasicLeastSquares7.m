clc
clear all

grav = 9.81;

ks__ = 500.4;
Bs__ = 24.67; 
Ms__ = 325;               % 1/4 sprung mass (kg)
Mus__ = 65;               % 1/4 unsprung mass (kg)
kus__ = 500.5e2;          % tire stiffness (N/m)
Bus__ = 0;

select = ([15:17 24 26]) % GOOD FOR ALL VALUES
% select = ([15:17 24]); Good for Ms_ Mus_ Bs_ kus_
%select = ([15:17 19 23:24]);
%select = ([15:20 23:24]);

%NOTE: Noise&Error is more pronounced when T is small
% how nosie is distributed makes a difference too, noise with non-random
% magintude is much easier


%% CONSTANTS
%T=0.01; %1 second between each 
%points = 2000; % MORE DATA POINTS RESULTS IN MORE ERROR FOR SIMULINK SYSTEM
T=0.1; %1 second between each 
points = 1000; % T
NOISE = 0.00;  

FILTER = 20;
len = points;
Simulation_Time = (points)*T;
time = (0:T:T*(points))';


%% MODELS


[A,B,C,D] = returnStateSpace(ks__, kus__, Bs__, Ms__, Mus__, Bus__)
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);

INPUTS = size(B,2);
OUTPUTS = size(C,2);
STATEVAR = size(A,1);

syms  Ms_ Mus_ Bs_ ks_  kus_
TAR = [Ms_ Mus_ Bs_ kus_ ks_];

[A_,B_,C_,D_] = returnStateSpace(ks_, kus_, Bs_, Ms_, Mus_, Bus__);
%Apply forward Euler Approximatio
[Aq_, Bq_,Cq_,Dq_] = FORWARDEULER(A_,B_,C_,D_,T);

%Check system is stable

eigA = eig(A);
for i=1:length(eigA)
    if eigA(i)>0
        EigenValues_Positive
    end
end



%% SIMULATE SYSTEM USING SIMULINK
% simin1 = [time u1];
% simin2 = [time u2];
out = sim('SIMPLE_RUN_MODEL_REAL2_U3');

X_C = [];
Y_C = [];
X_D = [];
Y_D = [];
Y_CP =[];
Y_CFP =[];
Y_CF = [];
X_CP_ = [];
U = [];

for i=1:INPUTS
    U = [U u3.Data(:,i)];
end

for i=1:STATEVAR
    %Contin Model
    Y_C = [Y_C CONTIN_Y.Data(:,i)];
    Y_CF = [Y_CF CONTIN_YF.Data(:,i)];
    
    %Actual Plant
    Y_CP = [Y_CP CONTIN_YP.Data(:,i)];
    Y_CFP = [Y_CFP CONTIN_YFP.Data(:,i)];
    X_CP_ = [ X_CP_ CONTIN_XP.Data(:,i)];
    
    % Discrete model
    Y_D = [Y_D DISCRETE_Y.Data(:,i)];
end


 
X_D =(inv(Cq_)*(Y_D'-Dq_*U'))';

X_C =(inv(C_)*(Y_C'-D_*U'))';
X_CF = (inv(C_)*(Y_CF'-D_*U'))';
X_CP = (inv(C_)*(Y_CP'-D_*U'))';
X_CFP = (inv(C_)*(Y_CFP'-D_*U'))';


%X_C =(inv(C)*(Y_C-D*U'))'
%% SOLVE USING SIMULATED SIMULINK SYSTEM


% Add noise variables based upon the inputs
%
syms N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14
N = [N1;N2;N3;N4;N5;N6;N7;N8;N9;N10;N11;N12;N13;N14];

diff = length(select)-length(TAR);
N = N(1:diff); %Error thrown when diff>length(N) need more symbols

X_data = X_CFP;
Y_data = Y_CFP;
[var_, val_, PSI] = LeastSquaresRegression3(X_data,U,select,TAR,Aq_,Bq_,N);

var_ = var_;

sol = solve(val_==var_,{TAR,N1});
ks_ = eval(sol.ks_);
kus_ = eval(sol.kus_);
Bs_ = eval(sol.Bs_);
Ms_ = eval(sol.Ms_);
Mus_ = eval(sol.Mus_);

%noise = [vpa(sol.N1) vpa(sol.N2) vpa(sol.N3) vpa(sol.N4) vpa(sol.N5) vpa(sol.N6)]



fprintf("SIMULINK SYSTEM\n")
fprintf("ks=%.3f ~ ks_=%.3f ---> Error:%.3f\n",ks__,ks_,((-(ks__-ks_)/ks__))*100)
fprintf("kus=%.3f ~ kus_=%.3f ---> Error:%.3f\n",kus__,kus_,((-(kus__-kus_)/kus__))*100)
fprintf("Bs=%.3f ~ Bs_=%.3f ---> Error:%.3f\n",Bs__,Bs_,((-(Bs__-Bs_)/Bs__))*100)
%fprintf("Bus=%.3f ~ Bus_=%.3f ---> Error:%.3f\n",Bus__,Bus_,((-(Bus__-Bus_)/Bus__))*100)
fprintf("Ms=%.3f ~ Ms_=%.3f ---> Error:%.3f\n",Ms__,Ms_,((-(Ms__-Ms_)/Ms__))*100)
fprintf("Mus=%.3f ~ Mus_=%.3f ---> Error:%.3f\n",Mus__,Mus_,((-(Mus__-Mus_)/Mus__))*100)

% Similarity between model & actual 
% figure
% plot(DISCRETE_Y.time,X_CP,DISCRETE_Y.time,X_C);
% title("Initial MOdel Vs Plant");
% legend(["x1","x2","x3","x4","x1_m","x2_m","x3_m","x4_m"])

% X is correctly extracted from Y
% figure
% plot(DISCRETE_Y.time,X_CP,DISCRETE_Y.time,X_CP_);
% title("Ensure Correct X");
% legend(["Actual X","Derived X"]);


figure
plot(DISCRETE_Y.time,X_CP,DISCRETE_Y.time,X_CFP);
title("Noise Efffect");
legend(["x1","x2","x3","x4","x1_f","x2_f","x3_f","x4_f"])



figure
plot(DISCRETE_Y.time,Y_data(:,[1 3]),DISCRETE_Y.time,U(:,1:2));
title("'continous simulink Y-PLANT");
legend(["y1","y3","Fc","zr"]);

figure
plot(DISCRETE_Y.time,X_data(:,1:2),DISCRETE_Y.time,U(:,1:2));
title("'continous simulink X -PLANT");
legend({'x1',"x2","Fc","Zr"});
%plot(NOISE_DATA.Time,NOISE_DATA.Data)
%legend({'Simulink Noise'})


%% Evaluate

% Change variables to check
%Ms_ = 369;
%Mus_ = 369;


[A,B,C,D] = returnStateSpace(ks_, kus_, Bs_, Ms_, Mus_, Bus__)
%Apply forward Euler Approximatio
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);


out2 = sim('SIMPLE_RUN_MODEL_REAL2_U3');

X_C2 = [];
Y_C2 = [];
X_D2 = [];
Y_D2 = [];
U2 = [];

for i=1:INPUTS
    U2 = [U2 u3.Data(:,i)];
end

for i=1:STATEVAR
    %X_C2 = [X_C2 CONTIN_X.Data(:,i)];
    Y_C2 = [Y_C2 CONTIN_Y.Data(:,i)];
end

for i=1:STATEVAR
    %X_D2 = [X_D2 DISCRETE_X.Data(:,i)];
    Y_D2 = [Y_D2 DISCRETE_Y.Data(:,i)];
end


X_D2 =(inv(Cq_)*(Y_D2'-Dq_*U2'))';
X_C2 =(inv(C_)*(Y_C2'-D_*U2'))';

figure
subplot(2,2,1);
plot(CONTIN_Y.time,X_data(:,1),DISCRETE_Y.time,X_D2(:,1));
title("X1")
legend({'Actual Continous System - IN1','Model Continous System - IN2','Discrete Model- IN2'});
subplot(2,2,2);
plot(CONTIN_Y.time,X_data(:,2),DISCRETE_Y.time,X_D2(:,2));
title("X2")
subplot(2,2,3);
plot(CONTIN_Y.time,X_data(:,3),DISCRETE_Y.time,X_D2(:,3));
title("X3")
subplot(2,2,4);
plot(CONTIN_Y.time,X_data(:,4),DISCRETE_Y.time,X_D2(:,4));
title("X4")

figure
plot(DISCRETE_Y.time,Y_C(:,[1 3]));
title("'continous simulink Y-Model");
legend(["y1","y3"]);

figure
plot(DISCRETE_Y.time,Y_CP(:,[1 3]));
title("'continous simulink Y-PLANT");
legend(["y1","y3"]);


error = 0;

minlen = min(length(X_C),length(X_D2))
for i=1:minlen 
    for j=1:size(X_C,2)
        error = error + abs(X_C(i,j) - X_D2(i,j));  
    end
end


fprintf("Total Error:%.3f\n",error)