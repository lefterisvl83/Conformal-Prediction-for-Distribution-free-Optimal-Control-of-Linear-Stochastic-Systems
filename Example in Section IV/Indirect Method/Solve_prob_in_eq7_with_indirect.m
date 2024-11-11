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
% disturbance parameters
M = 100; % number of samples
mu1 = -.01;
sigma1 = sqrt(0.005); % what is inside sqrt is meant to be variance
shape_g= 5.5;
theta_g=0.005;
%% design parameters obtained
% inidrect method
Kind = [-1.4140 -2.3412]; % state-feedback gain obtained from solving (22)
Phi_ind=[3.4644, 3.8069;3.8069, 5.6494]; % matrix shaping the RPI obtained from solving (22)
EKbound_ind=0.4408; % Obtained by Lemma 3 using Kind obtained by solving (22)
%
%% Solve (7) by tighting the constraints via direct method
%% constraints X_t
% Generate points on a unit circle
theta_gwnia = linspace(0, 2*pi, 100);
unit_circle = [cos(theta_gwnia); sin(theta_gwnia)]; % Points on the unit circle
% Constraint-Tighetening (Direct Method)
Xt = Polyhedron(xball_rad*unit_circle');
Et = Polyhedron(draw_ellipsis(Phi_ind)');
Zt = Xt-Et;
% variables
z = sdpvar(2,N+1,'full');
v = sdpvar(N, 1,'full'); 
% constraints and tightening
F = [z(:,1)==x0];
for i=2:N+1
    F = [F, v(i-1)>=-Q+EKbound_ind, v(i-1)<=Q-EKbound_ind];
    F = [F, z(:,i)==A*z(:,i-1)+B*v(i-1)];
    F = [F, Zt.A*z(:,i)<=Zt.b];
end
% objective
QQ=eye(N);
for i=1:N
    QQ(i,i)=1;
end
obj = 1*v'*QQ*v+0*reshape(z,[2*(N+1) 1])'*reshape(z,[2*(N+1) 1])+100*z(:,N+1)'*z(:,N+1);
% optimize
ops = sdpsettings('solver', 'mosek');
result=optimize(F,obj,ops);
zz=value(z);
vv=value(v);

%% disturbance samples
M=10000;
w1traj_samples = [];
w2traj_samples = [];
%
for j=1:M
    w1traj_samples = random('Normal', mu1, sigma1, 1, N);
    w2traj_samples = gamrnd(shape_g, theta_g, 1, N).*(randi([0, 1], 1, N)*2-1);
    wtraj_samples{j}=[w1traj_samples;w2traj_samples];
end
%
%
xx=cell(1,M);
uu=cell(1,M);
for j=1:M
    xx{j}(:,1)=x0;
    for i=1:N
        uu{j}(i)=Kind*(xx{j}(:,i)-zz(:,i))+vv(i);
        xx{j}(:,i+1)=A*xx{j}(:,i)+B*uu{j}(i)+wtraj_samples{j}(:,i);
    end
end
%% plots
figure;
hold on;
for j=1:M
    plot(xx{j}(1,:), xx{j}(2,:), '-', 'LineWidth', 2, 'Color', [1, 0.5, 0.5, .2]);
end
%% plot deterministic trajectory
plot(zz(1,:),zz(2,:),'-','LineWidth',1,'Color',[1, 0, 0, 1])
Xt.plot('alpha',0.1,'color','gray')
Zt.plot('alpha',0.,'color','green', 'LineStyle', '--')
Et.plot('alpha',0.,'color','blue', 'LineStyle', '--')
xlabel('$x_1$','FontSize',18,'Interpreter','Latex')
ylabel('$x_2$','FontSize',18,'Interpreter','Latex')
%% count probabilities
count_e=0;
count_u=0;
for j=1:M
    count_ei=0;
    for i=2:N+1
        if xx{j}(:,i)'*Pt*xx{j}(:,i)<=1 % Xt.A*xx{j}(:,i)<=Xt.b
            count_ei=count_ei+1;
        end
    end
    if count_ei==N
        count_e=count_e+1;
    end
end
for j=1:M
    count_ui=0;
    for i=1:N
        if Q*uu{j}(i)*uu{j}(i)<=1
            count_ui=count_ui+1;
        end
    end
    if count_ui==N
        count_u=count_u+1;
    end
end
count_e=count_e/M;
count_u=count_u/M;
result
[count_e count_u]


function ellipsoid_points=draw_ellipsis(Y)
% draw the ellipsoid w'Yw<1
L = chol(Y, 'lower'); 

% Generate points on a unit circle
theta_gwnia = linspace(0, 2*pi, 100);
unit_circle = [cos(theta_gwnia); sin(theta_gwnia)]; % Points on the unit circle

% Map the unit circle to the ellipsoid using the inverse of L
ellipsoid_points = L \ unit_circle;
end
