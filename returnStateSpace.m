% Return Symbolic or Valued Variables
function [A,B,C,D] = returnStateSpace(ks, kus, Bs, Ms, Mus, Bus)
A = [0 0 1 0;
    0 0 0 1;
    (-ks-kus)/Mus ks/Mus (-Bs-Bus)/Mus Bs/Mus;
    ks/Ms -ks/Ms Bs/Ms -Bs/Ms];
B = [ 0 0 0;
    0 0 0;
    1/Mus kus/Mus Bus/Mus
    1/Ms 0 0];
C = [-1 1 0 0;
    0 0 0 1;
    1 0 0 0;
    0 0 1 0];
D = [0 0 0;
    0 0 0;
    0 -1 0;
    0 0 0];
end