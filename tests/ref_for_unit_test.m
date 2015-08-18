Q1_tilde = [7/27, -4/27; -4/27, 10/27];
% Q1_tilde = {7/27, -4/27, -4/27, 10/27};

B = [-3; 1];
sol_H = 2\B';
sol_H2 = 2\sol_H
sol = B*sol_H2 + Q1_tilde