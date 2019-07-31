% Compute the essential spectrum using continuation methods outlined in
% Rademacher, Scheel, and Sandstede 2007
% Modify continuation (number of steps, branch, etc) through contPar struc

close all; clear all;

%% Select system
% Modify for Rossler or Karma
file_names.wt = '../data_files/Rossler_1D_wave_train_c2.mat';  % Asymptotic wave train
file_names.out_name = 'essential_spectrum_data/rossler_essSpec_c2.mat';

file_names.ess_spec_solver = 'Rossler_essSpec_fcn';    % Solve essential spectrum function: spiral dispersion relation
file_names.problem_solver_1D = 'Rossler_1D';            % Solve the 1D system (used for initial point)

%% Set up
wt = load(file_names.wt);
par = wt.par;
numPar = wt.numPar;
u_infty = wt.U;         % Asymptotic wave train

% Continuation parameters for essential spectrum continuation -- This is a SIMPLE
% continuation!
contPar.ds   = 0.1*1i;                % continuation step size: sign determines direction - should be imaginary since nu = i*gamma is treated as a complex number
contPar.Name = 'nu';   
contPar.free = 'lambda';
contPar.numContSteps = 50;     % Number of continuation steps (in each direction)
contPar.full_branch = 1;        % 0: only continues in direction of contPar.ds. 1: continues in both directions to obtain the full branch
contPar.branch = 0;             % Which branch to compute: 0 at origin, 2 +: based on order of eigenvalues


% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...  % Options for fsolve
                    'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500);


addpath ../Util/

ess_spec_solver_fh = str2func(file_names.ess_spec_solver);
fh_1D = str2func(file_names.problem_solver_1D);

%% Find initial point for continuation

% Compute linear operator
[L1, L2] = ComputeLinearOperator_1D(par,numPar); % 1st and 2nd derivative matrices in 1D

if contPar.branch > 0
    phase_cond.u_old = u_infty(1:numPar.nx);
    [f,J] = fh_1D([u_infty;par.omega],L1,L2,par,numPar,phase_cond);

    J = full(J(1:3*numPar.nx,1:3*numPar.nx));
    [vecs,val] = eig(J);
    val = diag(val);
    [val_sort, index] = sort(abs(val),'ascend');


    par.lambda = val(index(contPar.branch));
    disp(par.lambda)
    U0 = vecs(1:par.numVars*numPar.nx,index(contPar.branch));


else
% %%% Derivative of wave train 

    U0 = zeros(par.numVars*numPar.nx,1);
    for j = 1:par.numVars
        tmp = L1*u_infty((j-1)*numPar.nx + 1:j*numPar.nx);
        tmp = tmp./norm(tmp);
        U0((j-1)*numPar.nx + 1:j*numPar.nx,1) = tmp;
        
    end
   
    par.lambda = 0;

end

par.nu = 0;
    
init_lambda = par.lambda;
init_nu = par.nu;

initial_sol = [U0; par.(contPar.free)]; % Initial eigenvectors
phase_cond.u_old = U0;

% Make sure eigenvector is correct

uout = fsolve(@(y) ess_spec_solver_fh(y,u_infty,L1,L2,par,numPar,phase_cond),initial_sol,options);
par.(contPar.free) = uout(end);
 
phase_cond.u_old = uout(1:par.numVars*numPar.nx);
initial_sol = uout;
%% Continuation

free_vec = zeros(contPar.numContSteps+1,1);
cont_vec = zeros(contPar.numContSteps+1,1);

free_vec(1) = par.lambda;
cont_vec(1) = par.nu;

for k = 1:contPar.numContSteps
    
    par.(contPar.Name) = par.(contPar.Name) + contPar.ds;   % Add step to continuation parameter
   
    uout = fsolve(@(y) ess_spec_solver_fh(y,u_infty,L1,L2,par,numPar,phase_cond),uout,options);
    par.(contPar.free) = uout(end);
    phase_cond.u_old = uout(1:par.numVars*numPar.nx);

    free_vec(k+1)   = par.(contPar.free);
    cont_vec(k+1)   = par.(contPar.Name);
    

end

%% Continuation in other direction to get full branch

if contPar.full_branch == 1
    free_vec2 = zeros(contPar.numContSteps,1);
    cont_vec2 = zeros(contPar.numContSteps,1);

    par.lambda = init_lambda;  % Reset to start at origin again
    par.nu = init_nu;  

    uout = initial_sol;
    phase_cond.u_old = U0;

    contPar.ds = -1*(contPar.ds);

    for k = 1:contPar.numContSteps
       
        par.(contPar.Name) = par.(contPar.Name) + contPar.ds;   % Add step to continuation parameter
   
        uout = fsolve(@(y) ess_spec_solver_fh(y,u_infty,L1,L2,par,numPar,phase_cond),uout,options);
    	par.(contPar.free) = uout(end);
        phase_cond.u_old = uout(1:par.numVars*numPar.nx);
     
        free_vec2(k)   = par.(contPar.free);
        cont_vec2(k)   = par.(contPar.Name);

    end

     
    ess_spec = [flip(free_vec); free_vec2];  % form full essential spectrum branch
    nu_vec   = [flip(cont_vec); cont_vec2];
     
     
 else
     ess_spec = free_vec;
     nu_vec   = cont_vec;
     
 end
%% Plot and save data


figure(1);   % Plot the essential spectrum branch
plot(ess_spec,'b-','linewidth',3); hold on;
ylabel('Im(\lambda)');  xlabel('Re(\lambda)'); 
title(['Essential Spectrum, Branch = ' num2str(contPar.branch)]); 
plot([0,0], ylim,'color','k','linewidth',2);
plot(xlim, [0,0],'color','k','linewidth',2);
box on;
set(gca,'fontsize',20,'linewidth',2);


%save(file_names.out_name,'ess_spec','nu_vec','par','numPar','contPar')


