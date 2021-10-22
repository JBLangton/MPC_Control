function [Aq, Bq,Cq,Dq] = Plant_Model_F_Euler(T,A,B,C,D)
[Aq, Bq,Cq,Dq] = FORWARDEULER(A,B,C,D,T);
end