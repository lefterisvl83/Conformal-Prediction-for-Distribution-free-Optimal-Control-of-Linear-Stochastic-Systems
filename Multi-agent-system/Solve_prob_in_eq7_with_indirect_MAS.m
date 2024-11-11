clear all
clc
close all
%
%
%% system parameters
rng(21) % rng(6,9,10) = 99% rng = 99% the best so far | rng(21)=97% best so 
A = kron(eye(6), [1 0.5;0 1]);
B = kron(eye(6), [0; 0.5]);
N = 100;
x0=kron(ones(6,1),[.5 -.5]');
xball_rad=sqrt(10);
Pt=kron(eye(6),eye(2)/(xball_rad^2));
Q = kron(eye(6),1);
% disturbance parameters
M = 100; % number of samples
mu1 = -.01;
sigma1 = sqrt(0.005); % what is inside sqrt is meant to be variance
shape_g= 5.5;
theta_g=0.005;
%% design parameters obtained
% inidrect method
Kind = kron(eye(6),[-1.4140 -2.3412]); % state-feedback gain obtained from solving (22)
Phi_ind=kron(eye(6),[3.4644, 3.8069;3.8069, 5.6494]); % matrix shaping the RPI obtained from solving (22)
EKbound_ind=0.4408; % Obtained by Lemma 3 using Kind obtained by solving (22)
%
%% Solve (7) by tighting the constraints via direct method
% variables
z = sdpvar(12,N+1,'full');
v = sdpvar(6, N,'full'); 
% constraints and tightening
F = [z(:,1)==x0];
for i=2:N+1
    F = [F, v(:,i-1)'*v(:,i-1)<=(1-EKbound_ind)^2];
%     F = [F, v(i-1)>=-Q+EKbound_ind, v(i-1)<=Q-EKbound_ind];
    F = [F, z(:,i)==A*z(:,i-1)+B*v(:,i-1)];
    F = [F, z(:,i)'*z(:,i)<= (1/sqrt(max(eig(Pt)))-1/sqrt(min(eig(Phi_ind))))^2];
end
% objective
QQ=eye(N);
for i=1:N
    QQ(i,i)=1;
end
obj = 1*reshape(v,[6*N 1])'*reshape(v,[6*N 1])+100*z(:,N+1)'*z(:,N+1);
% optimize
ops = sdpsettings('solver', 'mosek');
tic
result=optimize(F,obj,ops);
toc
zz=value(z);
vv=value(v);
%
result
