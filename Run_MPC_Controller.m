
% Sample Time
T = 0.05;
NOISE = 0.05;
sim_time = 10;
N = 10;


beta = 10;
Qa = 5*beta; % for x2d
Qb = 0*beta; % for
Qc = 0*beta;
alpha = 10; % for R

%Calculate MPC
CalculateMPC

t = zeros(1,sim_time/T);
U = zeros(1,sim_time/T);
t(1) = 0;
U(1) = 0;

freq = 1/5;
Amp = 0.1;

for i=2:sim_time/T +1
    t(i) = (i-1)*T;
    U(i) = Amp*sin(2*pi*freq*t(i));  
end
plot(t,U)

INPUTS = timeseries(U,t)
sim('MPC_Controller_Sim')

h = figure
plot(Y.Time,Y.Data)
title('Output Over Time')
xlabel('Times(s)')
ylabel('Time(s)')
legend({'x2 - x1', 'dot(x2)',  'x1-zr',  'dot(x1)', 'Z_r'})

x1_error = sum(abs(Y.Data(:,1)),"All");
x2_error = sum(abs(Y.Data(:,2)),"All");
x3_error = sum(abs(Y.Data(:,3)),"All");

fprintf("\nError for x_1 x_2 x_3 = %.3f %3.f %3.f\n",x1_error,x2_error,x3_error)

%filename = sprintf('Fig_F%2f_A%2f_%d_%d_%d_%d'+"type"+".png", freq,Amp,Qa,Qb,Qc,alpha)
%saveas(h,filename)

