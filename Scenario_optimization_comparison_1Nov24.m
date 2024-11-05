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
M = 1000;% ceil((2/(1-theta))*((N^2-1)*log(2)-log(1e-3))); % sample size according to [9] eq. (8)
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
    j
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

figure;
hold on;
for j=1:M
    plot(value(x(1,:,j)), value(x(2,:,j)), '-', 'LineWidth', 2, 'Color', [1, 0.5, 0.5, .2]);
end

% Display memory information
whos
memInfo = memory;
fprintf('Total Memory: %.2f MB\n', memInfo.MemAvailableAllArrays / (1024^2));
fprintf('Used Memory: %.2f MB\n', (memInfo.MemUsedMATLAB) / (1024^2));

% %% disturbance samples
% w1traj_samples = [];
% w2traj_samples = [];
% %
% for j=1:M
%     w1traj_samples = random('Normal', mu1, sigma1, 1, N);
%     w2traj_samples = gamrnd(shape_g, theta_g, 1, N).*(randi([0, 1], 1, N)*2-1);
%     wtraj_samples{j}=[w1traj_samples;w2traj_samples];
% end
% %
% %
% xx=cell(1,M);
% uu=cell(1,M);
% for j=1:M
%     xx{j}(:,1)=x0;
%     for i=1:N
%         uu{j}(i)=Kind*(xx{j}(:,i)-zz(:,i))+vv(i);
%         xx{j}(:,i+1)=A*xx{j}(:,i)+B*uu{j}(i)+wtraj_samples{j}(:,i);
%     end
% end
% %% plots
% 
% %% plot deterministic trajectory
% plot(zz(1,:),zz(2,:),'-','LineWidth',1,'Color',[1, 0, 0, 1])
% % scatter(zz(1,:),zz(2,:),20, 'black', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
% Xt.plot('alpha',0.1,'color','gray')
% Zt.plot('alpha',0.,'color','green', 'LineStyle', '--')
% % Xt1_20.plot('alpha',0.1,'color','green','edgealpha',0.)
% % Xt21_40.plot('alpha',0.1,'color','green','edgealpha',0.)
% % Xt31_100.plot('alpha',0.1,'color','green','edgealpha',0.)
% % Eet1_20.plot('alpha',0.0,'color',[0.5 0.5 0.5], 'LineStyle', '--')
% % Eet21_40.plot('alpha',0.0,'color',[0.5 0.5 0.5], 'LineStyle', '--')
% % Eet31_100.plot('alpha',0.0,'color',[0.5 0.5 0.5], 'LineStyle', '--')
% xlabel('$x_1$','FontSize',18,'Interpreter','Latex')
% ylabel('$x_2$','FontSize',18,'Interpreter','Latex')
% %% count probabilities
% count_e=0;
% count_u=0;
% for j=1:M
%     count_ei=0;
%     for i=2:N+1
%         if xx{j}(:,i)'*Pt*xx{j}(:,i)<=1 % Xt.A*xx{j}(:,i)<=Xt.b
%             count_ei=count_ei+1;
%         end
%     end
%     if count_ei==N
%         count_e=count_e+1;
%     end
% end
% for j=1:M
%     count_ui=0;
%     for i=1:N
%         if Q*uu{j}(i)*uu{j}(i)<=1
%             count_ui=count_ui+1;
%         end
%     end
%     if count_ui==N
%         count_u=count_u+1;
%     end
% end
% count_e=count_e/M;
% count_u=count_u/M;
% [count_e count_u]


function ellipsoid_points=draw_ellipsis(Y)
% draw the ellipsoid w'Yw<1
L = chol(Y, 'lower'); 

% Generate points on a unit circle
theta_gwnia = linspace(0, 2*pi, 100);
unit_circle = [cos(theta_gwnia); sin(theta_gwnia)]; % Points on the unit circle

% Map the unit circle to the ellipsoid using the inverse of L
ellipsoid_points = L \ unit_circle;
end
