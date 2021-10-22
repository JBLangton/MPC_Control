%function ParsePrediction (XPredict)

T = 0.05;
TIME = 7;
iteration = cast(TIME/T,'uint8');

X = [1 1 1 1 1];

x1 = XPredict.Data(iteration,1:6:end)
x2 = XPredict.Data(iteration,2:6:end)
x3 = XPredict.Data(iteration,3:6:end)
x4 = XPredict.Data(iteration,4:6:end)
x5 = XPredict.Data(iteration,5:6:end)
figure
plot(7+XPredict.Time(1:N,1)',x5)
