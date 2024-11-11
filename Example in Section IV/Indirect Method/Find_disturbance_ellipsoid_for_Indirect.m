clear all
clc
close all
%
%% TRAINING
% distubance samples
rng(1)
theta=0.95;
N = 100; % Number of time steps
M = 100; % number of training samples
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
% introduce variables
Yhat=sdpvar(2);
% constraints
F = [Yhat>=0];
for j=1:M
    for i=1:N
        F=[F, norm(Yhat*wtraj_samples{j}(:,i),2)<=1];
    end
end

optimize(F,-log(det(Yhat)),sdpsettings('solver','sedumi'))
Yel=value(Yhat)'*value(Yhat);
%
ellipsoid_points=draw_ellipsis(Yel);
% Plot the ellipsoid

%% CALIBRATE - CONFORMALIZE
% new calibration disturbance samples
rng(2)
w1traj_samples = [];
w2traj_samples = [];
%
for j=1:M
    w1traj_samples = random('Normal', mu1, sigma1, 1, N);
    w2traj_samples = gamrnd(shape_g, theta_g, 1, N).*(randi([0, 1], 1, N)*2-1);
    wtraj_samples{j}=[w1traj_samples;w2traj_samples];
end
%
rng(1)
for j=1:M
    ri=[];
    for i=1:N
        ri=[ri;sqrt(wtraj_samples{j}(:,i)'*Yel*wtraj_samples{j}(:,i))];
    end
    Rwi(j)=max(ri);
end
Cw=quantile(Rwi,theta);
%
Y=Yel/(Cw*Cw) % This leads to the ellipsoid: w'Y w < 1 with CP



%--------------------------------------------
function ellipsoid_points=draw_ellipsis(Y)
% draw the ellipsoid w'Yw<1
L = chol(Y, 'lower'); 

% Generate points on a unit circle
theta_gwnia = linspace(0, 2*pi, 100);
unit_circle = [cos(theta_gwnia); sin(theta_gwnia)]; % Points on the unit circle
[V, D] = eig(Y);

% Map the unit circle to the ellipsoid using the inverse of L
ellipsoid_points = V * sqrt(inv(D)) * unit_circle;
end


