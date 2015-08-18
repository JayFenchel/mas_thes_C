Phi_bl = [5 2 0; 2 4 -1; 0 -1 3];
C = chol(Phi_bl)';
help = C\eye(3);
Phi_bl_I = C'\help

f2 = Phi_bl_I(1:2, 1:2)

Phi_bl = [5 2 0; 2 4 -1; 0 -1 2];
C = chol(Phi_bl)';
help = C\eye(3);
Phi_bl_I = C'\help

Q1_tilde = [7/27, -4/27; -4/27, 10/27];
% Q1_tilde = {7/27, -4/27, -4/27, 10/27};

A = [1, 2; 1, 2];
B = [-3; 1];
mAB = -[A';B']

AB = [A,B]

Phi_bl_I*mAB

sol_H = 2\B';
sol_H2 = 2\sol_H;
sol = B*sol_H2 + Q1_tilde;

Phi_bl_I*AB'
AB*Phi_bl_I*AB'


f = [4 -3; -3 4];
fc = chol(f)';

help = fc\eye(2);
solu = fc'\help

AB*Phi_bl_I*AB'+f2