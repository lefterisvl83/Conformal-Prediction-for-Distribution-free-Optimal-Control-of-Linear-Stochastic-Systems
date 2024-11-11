clear all
clc
close all
%
%
%% system parameters
rng(21)
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
sigma1 = sqrt(0.005); 
shape_g= 5.5;
theta_g=0.005;
%% design parameters obtained from GA and Lemma 3
% direct method
Kdir = [-0.24075159 -0.78717577]; % state-feedback gain obtained from GA
Ebound_dir=0.5785; % obtained by Lemma 3
EKbound_dir=0.1271; % obtained by Lemma 3
%
%% Solve (7) by tighting the constraints via direct method
% variables
z = sdpvar(2,N+1,'full');
v = sdpvar(N, 1,'full'); 
% constraints and tightening
F = [z(:,1)==x0];
for i=2:N+1
    F = [F, v(i-1)>=-Q+EKbound_dir, v(i-1)<=Q-EKbound_dir];
    F = [F, z(:,i)==A*z(:,i-1)+B*v(i-1)];
    F = [F, z(:,i)'*z(:,i)<=(xball_rad-Ebound_dir)^2];
end
% objective
QQ=eye(N);
for i=1:N
    QQ(i,i)=1/1;
end
obj = 1*v'*QQ*v+0*reshape(z,[2*(N+1) 1])'*reshape(z,[2*(N+1) 1])+100*z(:,N+1)'*z(:,N+1);
% optimize
ops = sdpsettings('solver', 'mosek');
result=optimize(F,obj,ops);
zz=value(z);
vv=value(v);

%% constraints X_t
% Generate points on a unit circle
theta_gwnia = linspace(0, 2*pi, 100);
unit_circle = [cos(theta_gwnia); sin(theta_gwnia)]; % Points on the unit circle
% Constraint-Tighetening (Direct Method)
Xt = Polyhedron(xball_rad*unit_circle');
Zt = Polyhedron((xball_rad-Ebound_dir)*unit_circle');
Et = Polyhedron(Ebound_dir*unit_circle');
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
        uu{j}(i)=Kdir*(xx{j}(:,i)-zz(:,i))+vv(i);
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
Et.plot('alpha',0.,'color', 'b', 'LineStyle','--')
xlabel('$x_1$','FontSize',18,'Interpreter','Latex')
ylabel('$x_2$','FontSize',18,'Interpreter','Latex')
%% count probabilities
count_e=0;
count_u=0;
for j=1:M
    count_ei=0;
    for i=2:N+1
        if xx{j}(:,i)'*Pt*xx{j}(:,i)<=1 
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
[count_e count_u]
result




