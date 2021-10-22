function [Aq, Bq,Cq,Dq] = Plant_Model_ZOH(T,A,B,C,D)
 sys = ss(A,B,C,D);
 sys = c2d(sys,T,'zoh');
 Aq = sys.A;
 Bq = sys.B;
 Cq = sys.C;
 Dq = sys.D;
end