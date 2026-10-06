clc;
clear;
close all;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% available PSDs
% N1: Normal, N = 100
% N2: Normal, N = 100
% N3: Normal, N = 100
% U1: Uniform, N = 100
% U2: Uniform, N = 100
% U3: Uniform, N = 100
% LTB1: Unimodal Left-tailed, N = 100
% LTB2: Unimodal Left-tailed, N = 100
% LTB3: Unimodal Left-tailed, N = 100
% RTB1: Unimodal Right-tailed, N = 100
% RTB2: Unimodal Right-tailed, N = 100
% RTB3: Unimodal Right-tailed, N = 100
% SB: Symmetric Bimodal, N = 100
% LTBIM1: Left-tailed bimodal, N = 100
% LTBIM2: Left-tailed bimodal, N = 100
% LTBIM3: Left-tailed bimodal, N = 100
% RTBIM1: Right-tailed bimodal, N = 100
% RTBIM2: Right-tailed bimodal, N = 100
% RTBIM3: Right-tailed bimodal, N = 100
% CCGN2A: Coarse coffee, N = 2(A)
% CCGN2B: Coarse coffee, N = 2(B)
% CCGN3: Coarse coffee, N = 3
% CCGN5: Coarse coffee, N = 5
% CCGN10: Coarse coffee, N = 10
% CCGN30: Coarse coffee, N = 30
% FCGN2A: Fine coffee, N = 2(A)
% FCGN2B: Fine coffee, N = 2(B)
% FCGN3: Fine coffee, N = 3
% FCGN5: Fine coffee, N = 5
% FCGN10: Fine coffee, N = 10
% FCGN30: Fine coffee, N = 30
% Mono: Monodisperse, N = 1

load('PSDs.mat');
PSD_label = 'CCGN5';
model.qj = PSD.(PSD_label).qj;
model.wj = PSD.(PSD_label).wj;
model.Nclasses = length(model.qj);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference value of Theta; Theta = 0 corresponds to batch operation
model.Theta = 0;
% gamm denotes the dimensionless concentration gamma used in the manuscript
model.gammIn = 0;      % Dimensionless inlet liquid-phase concentration
model.BiRef = 1000;     % Reference Biot number
model.BiMode = 1;       % 1: constant Bi; 2: size-dependent Bi
model.K = 1;            % Equilibrium partition coefficient
model.alpha = 1;        % Solution-to-solid volume ratio
% Solution-to-solid ratio for each size class
model.aj = model.alpha./model.wj;    
model.ReExp = 1/3;      % Reynolds-number exponent in the Sh correlation

if model.BiMode == 1
    model.Bi = model.BiRef*ones(size(model.qj));
else
    model.Bi = model.qj.^model.ReExp * model.BiRef;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.Nnodes = 100;

% Spatial discretization
model.xi   = linspace(0,1,model.Nnodes);
model.xi2  = model.xi.^2;
model.dxi  = 1/(model.Nnodes-1);
model.dxi2 = model.dxi^2;

% Number of equations
if model.Theta == 0
    model.Neq = model.Nnodes*model.Nclasses+1;
else
    model.Neq = model.Nnodes*model.Nclasses+2;
end

% Diffusion coefficient factor
model.factor = 1./(model.qj.^2 * model.dxi2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.y0 = zeros(model.Neq,1);  % Initial condition
model.y0(1:(model.Nnodes*model.Nclasses)) = 1;
model.tauSpan = [0 0.6];        % Dimensionless simulation time span

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tau,psi,psim,psiS,gamm,eta] = PSD_model(model);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot solutions

% Plot overall average solid-phase concentration, semilog plot
figure;
semilogy(tau,psiS);
xlabel('\tau');
ylabel('\psi_S');
grid on;
title(sprintf('Overall solid-phase concentration: $\\alpha = %.2g$, $K = %.2g$, $PSD = %s$', ...
    model.alpha, model.K, PSD_label),'Interpreter','latex');

% Plot overall average solid-phase concentration, linear plot
figure;
title('Overall solid-phase concentration');
plot(tau,psiS);
xlabel('\tau');
ylabel('\psi_S');
grid on;
title(sprintf('Overall solid-phase concentration: $\\alpha = %.2g$, $K = %.2g$, $PSD = %s$', ...
    model.alpha, model.K, PSD_label),'Interpreter','latex');

% Plot extraction yield
figure;
plot(tau,eta);
xlabel('\tau');
ylabel('\eta');
grid on;
title(sprintf('Extraction yield: $\\alpha = %.2g$, $K = %.2g$, $PSD = %s$', ...
    model.alpha, model.K, PSD_label),'Interpreter','latex');

% Plot liquid-phase concentration
figure;
plot(tau,gamm);
xlabel('\tau');
ylabel('\gamma');
grid on;
title(sprintf('Liquid-phase concentration: $\\alpha = %.2g$, $K = %.2g$, $PSD = %s$', ...
    model.alpha, model.K, PSD_label),'Interpreter','latex');

% Plot average solutions for selected classes
id = [3 4]; % Maximum value should not exceed model.Nclasses
figure;
plot(tau,psim(:,id));
xlabel('\tau');
ylabel('$\bar{\psi}_j$','Interpreter','latex');
legend(compose('q_j = %.3f',model.qj(id)),'Location','best');
grid on;
title(sprintf('Volume-average solutions for selected classes: $\\alpha = %.2g$, $K = %.2g$, $PSD = %s$', ...
    model.alpha, model.K, PSD_label),'Interpreter','latex');

% Plot concentration profiles for desired classes
ti  = linspace(min(tau),max(tau),150);
xii = linspace(0,1,50);

[timei,posi] = meshgrid(ti,xii);

for j = 1:length(id)

    psii = interp2(model.xi,tau,psi(:,:,id(j)),posi,timei);

    figure;
    surf(timei,posi,psii,'EdgeColor','none','FaceAlpha',0.5);

    xlabel('$\tau$','Interpreter','latex');
    ylabel('$\xi_j$','Interpreter','latex');
    zlabel('$\psi_j$','Interpreter','latex');

    title(sprintf('Concentration profile for size class $j = %d$, $q_j = %.3f$', ...
        id(j), model.qj(id(j))), ...
        'Interpreter','latex');

end