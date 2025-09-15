clear
clc
close all



A = [1,2;2,3];
B = [1,2;-1,-4];

num_states = size(A,2);
num_controls = size(B,2);

C = eye(size(A));
D = zeros(size(B));

Q = eye(size(A));
Q = diag([1 10]);
R = eye(size(B));

K = lqr(A,B,Q,R);

% Steady state ( x-->x_ss , u-->y_ss , y --> y_ss = r_ss )
% 0 = Ax + Bu
% y = Cx + Du

% let x_ss = N_x * r_ss
% let u_ss = N_u * r_ss

% generally, u = u_ss - K(x-x_ss)

%  [ A , B ; C , D ] * [N_x ; N_u ] = [ 0 ; 1 ] 

% u = Nu * r - K ( x - N_x *r )

% alternatively

% u = -K * x + N_bar * r
% N_bar = (N_u + KN_x)

G = [A,B;C,D];


r = eye(size(C,1));
N = G \ [zeros(size(A,1), size(C,1)); r];

N_x = N(1:size(A,1), :);
N_u = N(size(A,1)+1:end, :);

N_bar = (N_u + K*N_x)

N_bar = -inv(C * inv(A - B*K) * B)










