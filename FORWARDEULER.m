function [Aq,Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T)
dimA = size(A);
dimA = dimA(1)

Aq = (eye(dimA)+A*T)
Bq = B*T
Cq = C;
Dq = D;
end