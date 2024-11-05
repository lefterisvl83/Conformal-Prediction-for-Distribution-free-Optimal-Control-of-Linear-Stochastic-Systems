clear all
close all
clc
%
% distubance samples
theta=0.95;
N = 100; % Number of time steps
M = 100; % number of samples
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
% define variables
Yhat=sdpvar(2);
F = [Yhat>=0];
for j=1:M
    for i=1:N
        F=[F, norm(Yhat*wtraj_samples{j}(:,i),2)<=1];
    end
end


