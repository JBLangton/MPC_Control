ks = 500.4;
bs = 24.67; 

Ms = 325;               % 1/4 sprung mass (kg)
Mus = 65;               % 1/4 unsprung mass (kg)
kus = 500.5e2;          % tire stiffness (N/m)
bus = 0;

A = [ 0 0 1 0;
    0 0 0 1;
    -(ks +kus)/Mus ks/Mus -(bs+bus)/Mus bs/Mus;
    ks/Ms -ks/Ms bs/Ms -bs/Ms];
B = [0 0 ;
    0 0;
    1/Mus kus/Mus;
    1/Ms 0];
C = [-1 1 0 0;
    0 0 0 1;
    1 0 0 0;
    0 0 1 0];
D =  [0 0;0 0; 0 -1; 0 0];

