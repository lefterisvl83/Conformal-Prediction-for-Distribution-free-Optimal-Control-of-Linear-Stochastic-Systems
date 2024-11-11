clear all
clc
close all
%
%
%% system parameters
rng(21) % rng(6,9,10) = 99% rng = 99% the best so far | rng(21)=97% best so 
A = [1 0.5;0 1];
B = [0; 0.5];
N = 100;
x0=[2 -1]';
xball_rad=sqrt(10);
Pt=eye(2)/(xball_rad^2);
Q = 1;
theta=0.95; % specification probability
% disturbance parameters
m=1;n=2;
d=m*n*(N-1)*N/2;
beta=0.6;
M = 100; %ceil((2/(1-theta))*((d-1)*log(2)-log(beta))); % minimum sample size according to [9] eq. (8)
mu1 = -.01;
sigma1 = sqrt(0.005); % what is inside sqrt is meant to be variance
shape_g= 5.5;
theta_g=0.005;
% disturbance samples
%% disturbance samples
w1traj_samples = [];
w2traj_samples = [];
%
for j=1:M
    w1traj_samples = random('Normal', mu1, sigma1, 1, N);
    w2traj_samples = gamrnd(shape_g, theta_g, 1, N).*(randi([0, 1], 1, N)*2-1);
    wtraj_samples{j}=[w1traj_samples;w2traj_samples];
end
%% Solve (3) by Scenario Optimization
% define variables
x = sdpvar(2,N+1,M,'full');
g = sdpvar(N,1,'full');
Kti = sdpvar(2,N,N,'full');
u = sdpvar(N, 1,M,'full'); 
% Esum = sdpvar(1);
% constraints
F = [];
for j=1:M
    j;
    F = [F, x(:,1,j)==x0];
    F = [F, u(1,1,j)==g(1)];
    for i=2:N+1
        F = [F, x(:,i,j)==A*x(:,i-1,j)+B*u(i-1,1,j)+wtraj_samples{j}(:,i-1)];
        F = [F, x(:,i,j)'*Pt*x(:,i,j)<=1];
    end
    for i=2:N
        ww = [];
        for k=2:i
            ww = [ww; wtraj_samples{j}(:,i-k+1)];
        end
        F = [F, u(i,1,j)==g(i)+reshape(Kti(:,i,2:i),[1 2*(i-1)])*ww];
        F = [F, Q*u(i,1,j)*u(i,1,j)<=1];
    end
end
% objective
QQ=eye(N);
Esum=[];
for j=1:M
    Esum=Esum+u(:,1,j)'*QQ*u(:,1,j)+100*x(:,N+1,j)'*x(:,N+1,j);
end
obj = Esum/M;
% optimize
ops = sdpsettings('solver', 'mosek');
result=optimize(F,obj,ops);
%
% new M for plots
M =100;
figure;
hold on;
for j=1:M
    plot(value(x(1,:,j)), value(x(2,:,j)), '-', 'LineWidth', 2, 'Color', [1, 0.5, 0.5, .2]);
end




