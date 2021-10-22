clc
clear all
Simulation_Time = 500; 
points = 10000; % T
T = Simulation_Time/points; % Read Data Points

ts = T;             % sampling time step  

% High fliter is no filter
FREQ = 99999; % low pass filter 100
FILTER  = FREQ;

NOISE = 0.05;
SEED =100;

% M

ks__ = 500.4;
Bs__ = 24.67; 

Ms__ = 325;               % 1/4 sprung mass (kg)
Mus__ = 65;               % 1/4 unsprung mass (kg)
kus__ = 500.5e2;          % tire stiffness (N/m)
Bus__ = 0;
grav = 9.81;            % acceleration of gravity (m/s^2)
%v = 10;                 % vehicle velocity (m/s)

                        
% Construct linear state space model 
Aqcar = [-Bs__/Ms__ -ks__/Ms__ Bs__/Ms__ ks__/Ms__;  1 0  0 0; ...
    Bs__/Mus__ ks__/Mus__ -(Bs__+Bus__)/Mus__ -(ks__+kus__)/Mus__; 0 0 1 0];
Bqcar = [1/Ms__ 0 -1/Mus__ 0;0 0 kus__/Mus__ 0]'; 
Cqcar = eye(4); 
Dqcar = zeros(2,4)'; 
grav = 9.81;  


% From converting measured outputs to state space variables within Plant
% Straight from known model
DqCAR =  [0 0;0 0; 0 -1; 0 0];
CqCAR = [-1 1 0 0;
    0 0 0 1;
    1 0 0 0;
    0 0 1 0];


out = sim('RLS_Plant');

%ouputs ofinterest
a11 = RP1.Data(end,1);
a12 = RP1.Data(end,2);
a13 = RP1.Data(end,3);
a14 = RP1.Data(end,4);
a15 = RP1.Data(end,5);
a16 = RP1.Data(end,6);

a21 = RP3.Data(end,1);
a22 = RP3.Data(end,2);
a23 = RP3.Data(end,3);
a24 = RP3.Data(end,4);
a25 = RP3.Data(end,5);
a26 = RP3.Data(end,6);

a31 = RP3.Data(end,1);
a32 = RP3.Data(end,2);
a33 = RP3.Data(end,3);
a34 = RP3.Data(end,4);
a35 = RP3.Data(end,5);
a36 = RP3.Data(end,6);

a41 = RP4.Data(end,1);
a42 = RP4.Data(end,2);
a43 = RP4.Data(end,3);
a44 = RP4.Data(end,4);
a45 = RP4.Data(end,5);
a46 = RP4.Data(end,6);

syms ks kus Mus Ms Bs
[A,B,C,D] = returnStateSpace(ks, kus, Bs, Ms, Mus, Bus__)
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);

INPUTS = size(B,2);
OUTPUTS = size(C,2);
STATEVAR = size(A,1);

var = [];
for i =1:size(Aq,1)
    var = [var ;transpose(Aq(i,:)) ;transpose(Bq(i,:))];
end

select = ([15:17 24 26])
var = var(select);% GOOD FOR ALL VALUES
% a31 a32 a33 a43 a45

val = [a31 a32 a33 a43 a45]';

sol = solve(val==var,[ks kus Mus Ms Bs])

ks = eval(sol.ks(1));
kus = eval(sol.kus(1));
Bs = eval(sol.Bs(1));
Ms = eval(sol.Ms(1));
Mus = eval(sol.Mus(1));


fprintf("SIMULINK SYSTEM: %d Solutions\n",length(sol.Ms))
fprintf("ks_=%.3f ~ ks=%.3f ---> Error:%.3f\n",ks__,ks,(((ks-ks__)/ks__))*100)
fprintf("kus_=%.3f ~ kus=%.3f ---> Error:%.3f\n",kus,kus__,(((kus-kus__)/kus__))*100)
fprintf("Bs_=%.3f ~ Bs=%.3f ---> Error:%.3f\n",Bs,Bs__,(((Bs-Bs__)/Bs__))*100)
fprintf("Ms_=%.3f ~ Ms=%.3f ---> Error:%.3f\n",Ms,Ms__,(((Ms-Ms__)/Ms__))*100)
fprintf("Mus_=%.3f ~ Mus=%.3f ---> Error:%.3f\n",Mus,Mus__,(((Mus-Mus__)/Mus__))*100)


param_error = ( abs(((ks-ks__)/ks__)*100) + abs(((kus-kus__)/kus__)*100) + abs(((Bs-Bs__)/Bs__)*100) ...
    + abs(((Ms-Ms__)/Ms__)*100) + abs(((Mus-Mus__)/Mus__))*100 )/5;

%% EVALUATE

%Actual Model
% Gettng Paramter and Subbing In
[A,B,C,D] = returnStateSpace(ks__, kus__, Bs__, Ms__, Mus__, Bus__)


%Apply forward Euler Approximatio
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);

% Getting State Space Directly
Aq = [RP1.Data(end,[1:4]);
    RP2.Data(end,[1:4]);
    RP3.Data(end,[1:4]);
    RP4.Data(end,[1:4])];

Bq = [RP1.Data(end,[5:7]);
    RP2.Data(end,[5:7]);
    RP3.Data(end,[5:7]);
    RP4.Data(end,[5:7])];


out2 = sim('simModel');


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


X_D2 =(inv(Cq)*(Y_D2'-Dq*U2'))';
X_C2 =(inv(C)*(Y_C2'-D*U2'))';

figure
subplot(2,2,1);
plot(CONTIN_Y.time,X_C2(:,1),DISCRETE_Y.time,X_D2(:,1));
title("X1")
legend({'Actual Continous System - IN1','Discrete Model- IN2'});
subplot(2,2,2);
plot(CONTIN_Y.time,X_C2(:,2),DISCRETE_Y.time,X_D2(:,2));
title("X2")
subplot(2,2,3);
plot(CONTIN_Y.time,X_C2(:,3),DISCRETE_Y.time,X_D2(:,3));
title("X3")
subplot(2,2,4);
plot(CONTIN_Y.time,X_C2(:,4),DISCRETE_Y.time,X_D2(:,4));
title("X4")


error = 0;

minlen = min(length(X_C2),length(X_D2))
loop = (minlen-1)*10/(Simulation_Time) +1
for i=1:loop
    for j=1:size(XP.Data,2)
        error = error + abs(X_C2(i,j) - X_D2(i,j));  
    end
end
error = error/loop;


fprintf("Sim Error:%.3f\n",error)
fprintf("Param Error:%.3f\n",param_error)


% See effect of filter
figure
plot(Output_P.Time,Output_P.Data,XP.Time,XP.Data)
legend({'x1f','x2f','x3f','x4f','x1','x2','x3','x4'})

A_discrete_model = Aq;
B_discrete_model = Bq;
C_discrete_model = Cq;
D_discrete_model = Dq;