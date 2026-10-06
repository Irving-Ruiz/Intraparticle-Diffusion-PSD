function [tau,psi,psim,psiS,gamm,eta] = PSD_model(model)

% Solve ODE system
options = odeset('RelTol',1E-6,'AbsTol',1E-9);
[tau,y] = ode15s(@(tau,y) RHS(tau,y,model), ...
    model.tauSpan,model.y0,options);

% Model dimensions
M  = model.Nnodes;
N  = model.Nclasses;
Nt = length(tau);

% Recover solid-phase solution
psi = y(:,1:M*N);
psi = reshape(psi,Nt,M,N);

% Recover liquid-phase variables
gamm = y(:,M*N+1);

% Average solid-phase concentration for each size class
u = psi.*model.xi2;

psim = zeros(Nt,N);

for j = 1:N
    psim(:,j) = 3*volave(u(:,:,j),model.dxi);
end

% Overall average solid-phase concentration
psiS = psim*model.wj;

% Extraction yield
if model.Theta == 0
    eta = (1-psiS)*model.K./(model.K+gamm);
else
    eta = y(:,M*N+2);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = RHS(tau,y,model)

% Model parameters
M        = model.Nnodes;
N        = model.Nclasses;
Neq      = model.Neq;
K        = model.K;
Bi       = model.Bi;
Theta    = model.Theta;
gammIn   = model.gammIn;
alpha    = model.alpha;
aj       = model.aj;
factor   = model.factor;
dxi      = model.dxi;
qj       = model.qj;

% Initialize RHS
f = zeros(Neq,1);

% Recover state variables
psi  = reshape(y(1:M*N),M,N);
gamm = y(M*N+1);
if Theta == 0
    gammb = K*psi(M,:)' + K/alpha;
else   
    gammb = K*psi(M,:)';
    % GammIn = 0 implies beta = K
    beta = K;   
end

% Diffusion equations
dpsi = zeros(M,N);

for j = 1:N   
    for k = 1:M  
        if k == 1            
            dpsi(k,j) = 2*psi(k+1,j) - 2*psi(k,j);
        elseif k == M           
            psiM1j = psi(k-1,j) - 2*dxi*Bi(j)*(gammb(j) - gamm)/K;
            dpsi(k,j) = psiM1j - 2*psi(k,j) + psi(k-1,j) ...
                + (psiM1j - psi(k-1,j))/(k-1); 
        else
            dpsi(k,j) = psi(k+1,j) - 2*psi(k,j) + psi(k-1,j) ...
                + (psi(k+1,j) - psi(k-1,j))/(k-1);
        end
        dpsi(k,j) = factor(j)*dpsi(k,j);
    end
end

% Liquid-phase balance
dgamm = sum(3*Bi.*(gammb-gamm)./aj./qj.^2)+Theta*(gammIn-gamm);
% Assemble RHS
if Theta == 0
    f(1:end-1) = dpsi(:);
    f(end) = dgamm;
else
    % Extraction yield equation
    deta  = alpha*Theta*(gamm-gammIn)/beta;

    f(1:end-2) = dpsi(:);
    f(end-1) = dgamm;
    f(end) = deta;
end

end

function up = volave(u,dxi)
up = sum((u(:,1:end-1)+u(:,2:end))*dxi/2,2);
end
