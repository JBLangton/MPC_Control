clc
clear all

Bus = 0;
ks = 100;
kus = 200000.1;

Bs = 2;
Ms = 55;
Mus = 100;

select = ([15:17 24 26]) % GOOD FOR ALL VALUES
% select = ([15:17 24]); Good for Ms_ Mus_ Bs_ kus_
%select = ([15:17 19 23:24]);
%select = ([15:20 23:24]);

%NOTE: Noise&Error is more pronounced when T is small
% how nosie is distributed makes a difference too, noise with non-random
% magintude is much easier


%% CONSTANTS
T=0.0001; %1 second between each 
points = 100000; % T
NOISE = 0.00  

len = points;
Simulation_Time = (points)*T; % The longer the simulation time the greater the error past a certain point
time = (0:T:T*(points))';
Simulation_Time2 = Simulation_Time; % The shorter the evaluated period the more accurate, as there are less internal dynamics that are modelled


%% MODELS


[A,B,C,D] = returnStateSpace(ks, kus, Bs, Ms, Mus, Bus)
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);

INPUTS = size(B,2);
OUTPUTS = size(C,2);
STATEVAR = size(A,1);

syms  Ms_ Mus_ Bs_ ks_  kus_
TAR = [Ms_ Mus_ Bs_ kus_ ks_];

[A_,B_,C_,D_] = returnStateSpace(ks_, kus_, Bs_, Ms_, Mus_, Bus);
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
out = sim('RUN_MODEL_U_3');

u3 = out.u3.Data;

X_C = [];
Y_CF = [];
Y_C = [];
X_D = [];
Y_D = [];
U = [];

for i=1:INPUTS
    U = [U out.u3.Data(:,i)];
end

for i=1:STATEVAR
    %X_C = [X_C out.CONTIN_X.Data(:,i)];
    Y_C = [Y_C out.CONTIN_Y.Data(:,i)];
    Y_CF = [Y_CF out.CONTIN_YF.Data(:,i)];
    Y_D = [Y_D out.DISCRETE_Y.Data(:,i)];
end

for i=1:STATEVAR
    %X_D = [X_D out.DISCRETE_X.Data(:,i)];
    
end

X_D =(inv(Cq_)*(Y_D'-Dq_*U'))';
X_C =(inv(C_)*(Y_C'-D_*U'))';
X_CF =(inv(C_)*(Y_CF'-D_*U'))';
%X_C =(inv(C)*(Y_C-D*U'))'
%% SOLVE USING SIMULATED SIMULINK SYSTEM

check = (Y_C'-D_*U')';
figure
plot(out.DISCRETE_Y.time,Y_C(:,3),out.DISCRETE_Y.time,U(:,2),out.DISCRETE_Y.time,check(:,3))
legend({'Y_C','U','Y_C-D_*U'});

% Add noise variables based upon the inputs
%
syms N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12 N13 N14
N = [N1;N2;N3;N4;N5;N6;N7;N8;N9;N10;N11;N12;N13;N14];

diff = length(select)-length(TAR);
N = N(1:diff); %Error thrown when diff>length(N) need more symbols


[var_, val_, PSI] = LeastSquaresRegression3(X_C ,U,select,TAR,Aq_,Bq_,N);

var_ = var_;% +N1;


sol = solve(val_==var_,TAR);
ks_ = eval(sol.ks_(1));
kus_ = eval(sol.kus_(1));
Bs_ = eval(sol.Bs_(1));
Ms_ = eval(sol.Ms_(1));
Mus_ = eval(sol.Mus_(1));

%noise = [vpa(sol.N1) vpa(sol.N2) vpa(sol.N3) vpa(sol.N4) vpa(sol.N5) vpa(sol.N6)]



fprintf("SIMULINK SYSTEM: %d Solutions\n",length(sol.Ms_))
fprintf("ks=%.3f ~ ks_=%.3f ---> Error:%.3f\n",ks,ks_,((-(ks-ks_)/ks))*100)
fprintf("kus=%.3f ~ kus_=%.3f ---> Error:%.3f\n",kus,kus_,((-(kus-kus_)/kus))*100)
fprintf("Bs=%.3f ~ Bs_=%.3f ---> Error:%.3f\n",Bs,Bs_,((-(Bs-Bs_)/Bs))*100)
%fprintf("Bus=%.3f ~ Bus_=%.3f ---> Error:%.3f\n",Bus,Bus_,((-(Bus-Bus_)/Bus))*100)
fprintf("Ms=%.3f ~ Ms_=%.3f ---> Error:%.3f\n",Ms,Ms_,((-(Ms-Ms_)/Ms))*100)
fprintf("Mus=%.3f ~ Mus_=%.3f ---> Error:%.3f\n",Mus,Mus_,((-(Mus-Mus_)/Mus))*100)

figure
plot(out.DISCRETE_Y.time,Y_C(:,1),out.DISCRETE_Y.time,Y_D(:,1));
legend({'continous simulink','discrete simulink'});

figure
plot(out.DISCRETE_Y.time,Y_C(:,[1 3]),out.DISCRETE_Y.time,U(:,1:2));
title("'continous simulink Y - MODEL");
legend({'y1',"y3","Fc",'Zr'});

figure
plot(out.DISCRETE_Y.time,X_C(:,[1 2]),out.DISCRETE_Y.time,U(:,1:2));
title("'continous simulink X - MODEL");
legend({'x1',"x2","Zr","Fc"});
% plot(out.NOISE.Time,out.NOISE.Data)
% legend({'Simulink Noise'})


%% Evaluate

% Change variables to check
%Ms_ = 369;
%Mus_ = 369;


[A,B,C,D] = returnStateSpace(ks_, kus_, Bs_, Ms_, Mus_, Bus)
%Apply forward Euler Approximatio
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);

Simulation_Time = Simulation_Time2;
out2 = sim('RUN_MODEL_U_3');

X_C2 = [];
Y_C2 = [];
X_D2 = [];
Y_D2 = [];
U2 = [];

for i=1:INPUTS
    U2 = [U2 out2.u3.Data(:,i)];
end

for i=1:STATEVAR
    %X_C2 = [X_C2 out2.CONTIN_X.Data(:,i)];
    Y_C2 = [Y_C2 out2.CONTIN_Y.Data(:,i)];
end

for i=1:STATEVAR
    %X_D2 = [X_D2 out2.DISCRETE_X.Data(:,i)];
    Y_D2 = [Y_D2 out2.DISCRETE_Y.Data(:,i)];
end


X_D2 =(inv(Cq_)*(Y_D2'-Dq_*U2'))';
X_C2 =(inv(C_)*(Y_C2'-D_*U2'))';


figure
subplot(2,2,1);
plot(out.CONTIN_Y.time,X_C(:,1),out2.CONTIN_Y.time,X_C2(:,1),out2.DISCRETE_Y.time,X_D2(:,1));
title("X1")
legend({'Actual Continous System - IN1','Model Continous System - IN2','Discrete Model- IN2'});
subplot(2,2,2);
plot(out.CONTIN_Y.time,X_C(:,2),out2.CONTIN_Y.time,X_C2(:,2),out2.DISCRETE_Y.time,X_D2(:,2));
title("X2")
subplot(2,2,3);
plot(out.CONTIN_Y.time,X_C(:,3),out2.CONTIN_Y.time,X_C2(:,3),out2.DISCRETE_Y.time,X_D2(:,3));
title("X3")
subplot(2,2,4);
plot(out.CONTIN_Y.time,X_C(:,4),out2.CONTIN_Y.time,X_C2(:,4),out2.DISCRETE_Y.time,X_D2(:,4));
title("X4")
%figure
%plot(out.CONTIN_Y.time,U(:,1),out2.DISCRETE_Y.time,U2(:,1));
%legend({'Input1 - Initital','Input1 - Test'});

error=0;
minlen = min(length(X_C),length(X_D2))
for i=1:minlen 
    for j=1:size(X_C,2)
        error = error + abs(X_C(i,j) - X_D2(i,j));  
    end
end


fprintf("Total Error:%.3f\n",error)

fprintf("Total Error:%.3f\n",error)