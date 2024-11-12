clear all
close all
clc
%
%% disturbance samples for Calibration in Lemma 3
rng(41)
theta=0.95;
N = 100; % Number of time steps
M = 100; %277508; % number of samples
mu1 = -.01;
sigma1 = sqrt(0.005); % what is inside sqrt is meant to be variance
shape_g= 5.5;
theta_g=0.005;

%
w1traj_samples = [];
w2traj_samples = [];
%

for j=1:M
    w1traj_samples = random('Normal', mu1, sigma1, 1, N);
    w2traj_samples = gamrnd(shape_g, theta_g, 1, N).*(randi([0, 1], 1, N)*2-1);
    wtraj_samples{j}=[w1traj_samples;w2traj_samples];
end
%
% system matrices
A = [1 0.5;0 1];
B = [0;0.5];
Kdir = [-0.24075159 -0.78717577]; % state-feedback gain obtained from GA
Kind = [-1.4140   -2.3412]; % state-feedback gain obtained from solving (22)
K = Kind; % select Kind
%
% error trajectory samples for calibration
for j=1:M
    re{j}=[];
    ru{j}=[];
    for i=1:N
        error_samples{j}(:,1)=zeros(2,1);
        error_samples{j}(:,i+1)=(A+B*K)*error_samples{j}(:,i)+wtraj_samples{j}(:,i);
        re{j}=[re{j} norm(error_samples{j}(:,i+1))];
        ru{j}=[ru{j} norm(K*error_samples{j}(:,i))];
    end
    Re(j)=max(re{j});
    Ru(j)=max(ru{j});
end
Ce=quantile(Re, theta);
Cu=quantile(Ru, theta);
[Ce Cu]
    








