Phi_bl = [5 2 0; 2 4 -1; 0 -1 2];
C = chol(Phi_bl)';
help = C\eye(3);
Phi_bl_I = C'\help

Q1_tilde = [7/27, -4/27; -4/27, 10/27];
% Q1_tilde = {7/27, -4/27, -4/27, 10/27};

A = [1, 2; 1, 2];
B = [-3; 1];
mAB = -[A';B']

Phi_bl_I*mAB

sol_H = 2\B';
sol_H2 = 2\sol_H;
sol = B*sol_H2 + Q1_tilde;