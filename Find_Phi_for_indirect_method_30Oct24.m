clear all
clc
close all
%
%% Find Kind
% system parameters
tic
A = [1 0.5;0 1];
B = [0; 0.5];
Pt = (1/9)*eye(2);
Y = inv([12.1932 -1.0314;-1.0314 110.4454]); % disturbance matrix
eps=1e-6;
Q = 1;
% solve (22)
step_size=0.05 ;
lambda0_vals = [];
lambda1_vals = [];
% Loop over possible values for lambda0_vals(i)
for lambda0 = 0.001:step_size:1
    % For each lambda0, calculate the possible range of lambda1_vals(i) (which must satisfy lambda0 + lambda1_vals(i) <= 1)
    for lambda1 = 0.001:step_size:(1 - lambda0)
        % Store the valid pair (lambda0, lambda1)
        lambda0_vals = [lambda0_vals; lambda0];
        lambda1_vals = [lambda1_vals; lambda1];
    end
end
MM=length(lambda1_vals);
Phihatvalue=cell(1,MM);
Psivalue=cell(1,MM);
for i=1:MM
    clear('Phihat', 'Psi');
    Phihat=sdpvar(2);
    Psi = sdpvar(1,2);
    F = [];
    F = [Phihat <= inv(Pt)+eye(2)*eps];
    F = [F, [Phihat-(1/lambda1_vals(i))*inv(Y), A*Phihat+B*Psi; (A*Phihat+B*Psi)', lambda0_vals(i)*Phihat]>=0];
    F = [F, [Phihat, Psi'*sqrt(Q); sqrt(Q)*Psi, eye(1)]>=0];
    F = [F, Phihat<=inv(Pt)+eps*eye(2)];
    F = [F, Phihat>=0];
    optimize(F,trace(Phihat));
    Phihatvalue{i} = value(Phihat);% [value(Phi), obj, lambda0_vals(i), lambda1_vals(i)];
    Psivalue{i} = value(Psi);
    obj_value(i)=trace(value(Phihat));
end
%
pickmatrix=[1 obj_value(1)];
for i=1:length(obj_value)
    if obj_value(i)>0 && obj_value(i)<pickmatrix(2)
        pickmatrix=[i obj_value(i)];
    end
end
Kind=Psivalue{pickmatrix(1)}*inv(Phihatvalue{pickmatrix(1)})
eig(Ak=A+B*Kind);
Phiopt=inv(Phihatvalue{pickmatrix(1)});
eig(Phiopt)
toc