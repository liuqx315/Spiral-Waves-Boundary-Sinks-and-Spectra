%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical continuation of spiral wave with natural boundary conditions
% on a bounded disk and associated eigenfunction
% On each step, parameter is updated, spiral solved, eigenfunction solved
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%% Select system to solve
% Modify for Rossler or Karma
file_names.spiral = '../data_files/Rossler_spiral_R125_c3p4.mat'; % Initial data
file_names.efcn = '../data_files/rossler_c3p4_lineDefect_efcn2.mat';  % Initial eigenfunction
file_names.out_name = '../data_files/outfile.mat';  % Outfile 

file_names.spiral_problem = 'Rossler_2D_spiral_nonreflecting';  % System to solve
file_names.efcn_problem = 'Rossler_nonReflectingBoundary_evalProb';  % System to solve


%% Set up
u0 = load(file_names.spiral);
U0 = u0.U;

par = u0.par;
numPar = u0.numPar;

efcn = load(file_names.efcn);
Ubar0 = efcn.Ubar;
par.lambda = efcn.par.lambda;

% Continuation parameter: picked from system parameter
contPar.final = 0.2014;
contPar.ds = 0.01;                  % continuation step size
contPar.Name = 'kappa';
contPar.Free1 = 'omega';    % Free variable for spiral solver
contPar.Free2 = 'lambda';   % Free variable for eigenfunction solver
contPar.numContSteps = 50;
contPar.plot_iter = 10;     % Continuation plots updated every plot_iter

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500);
      
addpath ../Util/

fh_spiral = str2func(file_names.spiral_problem); 
fh_eval   = str2func(file_names.efcn_problem); 

%% Set up differentiation matrices, boundary conditions, phase conditions

%Phase condition: picks out R/2
pc = zeros(numPar.nx,numPar.ny);
pc(:,ceil(numPar.ny/2)) = 1; 
phase_cond.pc = pc(:); %reshape(pc,numPar.nx*numPar.ny,1);

% Boundary condition matrices
[Ltheta,~] = ComputeLinearOperator_1D(par,numPar);  % 1D theta differentiation matrix

h = 1/(numPar.ny-1);
Dr = zeros(1,numPar.ny);
Dr(end-4:end) = [1/4; -4/3; 3; -4; 25/12]./h;
Ix = speye(numPar.nx,numPar.nx);
Lr = kron(Dr,Ix);

% Short grid
Dr = zeros(1,numPar.ny-1);
Dr(end-4:end) = [1/4; -4/3; 3; -4; 25/12]./h;
Ix = speye(numPar.nx,numPar.nx);
Lr_short = kron(Dr,Ix);
Lr_short = [zeros(numPar.nx,1), Lr_short];
Lr_short = sparse(Lr_short);

numPar.bndry_idx = (numPar.nx*(numPar.ny-1)+1):numPar.nx*numPar.ny;

[L1, L2] = ComputeLinearOperator(par,numPar);

numPar.L1 = L1; numPar.L2 = L2;
numPar.Ltheta = Ltheta; numPar.Lr = Lr;
%% Compute First point

% Spiral 
% Phase condition
us = (phase_cond.pc).*U0(1:numPar.nx*numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);
us_th = (L1*U0(1:numPar.nx*numPar.ny)); % d/d(theta) of phase condition
phase_cond.u_star_th = us_th(phase_cond.pc > 0);

phase_cond.Ubarold = Ubar0(1:numPar.nx*numPar.ny);

initial_sol = [U0; par.(contPar.Free1)];
initial_sol_efcn = [Ubar0;  par.(contPar.Free2)];

uout_sp = fsolve(@(y) fh_spiral(y,par,numPar, phase_cond),initial_sol,options); 
par.(contPar.Free1) = uout_sp(end);

uout_efcn = fsolve(@(y) fh_eval(y,uout_sp,par,numPar, phase_cond),initial_sol_efcn,options); 
par.(contPar.Free2) = uout_efcn(end);

%% Continuation data
free_vec1 = zeros(contPar.numContSteps,1);
free_vec2 = zeros(contPar.numContSteps,1);
cont_vec = zeros(contPar.numContSteps,1);
free_vec1(1) = par.(contPar.Free1);
free_vec2(1) = par.(contPar.Free2);
cont_vec(1) = par.(contPar.Name);

cont_efcn = zeros(3*numPar.nx*numPar.ny,contPar.numContSteps);
cont_efcn(:,1) = uout_efcn(1:3*numPar.nx*numPar.ny);

%% Continuation - Simple
con_start = tic;

for k = 1:contPar.numContSteps
    
    par.(contPar.Name) = par.(contPar.Name) + contPar.ds;   % Add step to continuation parameter

    disp(k)
    disp(['Continuation Parameter: ' num2str(par.(contPar.Name))]);
    
    % Solve spiral problem: Need new spiral for eigenfunction linearization
    uout_sp = fsolve(@(y) fh_spiral(y,par,numPar, phase_cond), uout_sp,options); 
    par.(contPar.Free1) = uout_sp(end); 
    
    % Solve eigenfunction problem
    uout_efcn = fsolve(@(y) fh_eval(y, uout_sp, par,numPar, phase_cond),uout_efcn,options); 
    par.(contPar.Free2) = uout_efcn(end); 
    
    % Update parameters
    free_vec1(k+1) = par.(contPar.Free1);
    free_vec2(k+1) = par.(contPar.Free2);
    cont_vec(k+1) = par.(contPar.Name);
    
    cont_efcn(:,k+1) = uout_efcn(1:3*numPar.nx*numPar.ny);
    
    % Prepare for next continuation step
    us = (phase_cond.pc).*uout_sp(1:numPar.nx*numPar.ny);  % phase condition for the next round of iteration
    phase_cond.u_star = us(phase_cond.pc>0);
    
    us_th = (L1*uout_sp(1:numPar.nx*numPar.ny));           % Will need to recompute linear operators if changing radius -- need to look into this (probably should be in the FixedPointSecantCorrector)
    phase_cond.u_star_th = us_th(phase_cond.pc > 0);
   
    phase_cond.Ubarold = uout_efcn(1:numPar.nx*numPar.ny);
   
    if mod(k,contPar.plot_iter) == 0 
         
        figure(1); plot(cont_vec(1:k+1),free_vec1(1:k+1),'-o','LineWidth',2); % continuation diagram
        xlabel(contPar.Name); ylabel(contPar.Free1); title('Continuation Parameter'); box on; set(gca,'fontsize',18,'linewidth',2); drawnow;
        
        figure(2); plot(cont_vec(1:k+1),real(free_vec2(1:k+1)),'-o','LineWidth',2); % continuation diagram
        xlabel(contPar.Name); ylabel(contPar.Free2); title('Continuation Parameter'); box on; set(gca,'fontsize',18,'linewidth',2); drawnow;
        
        
        % Spiral
        plot_spiral(uout_sp,par,numPar);
        
    end 
    
    % Stopping criteria: 
    if (sign(contPar.ds)*par.(contPar.Name)) >= (contPar.final*sign(contPar.ds))
        disp(['Continuation parameter reached before specified number of iterations. Final parameter value: ' num2str(par.(contPar.Name))])
            
        % Final continuation plot
        cont_vec = cont_vec(1:k+1);
        free_vec1 = free_vec1(1:k+1);
        free_vec2 = free_vec2(1:k+1);
        cont_efcn = cont_efcn(:,1:k+1);

        break;
    end
        
           
end

if k < contPar.numContSteps
    % One final solve to set parameter values
    par.(contPar.Name) = contPar.final;
    
    uout_sp = fsolve(@(y) fh_spiral(y,par,numPar, phase_cond),uout_sp,options); 
    par.(contPar.Free1) = uout_sp(end);

    uout_efcn = fsolve(@(y) fh_eval(y,uout_sp,par,numPar, phase_cond),uout_efcn,options); 
    par.(contPar.Free2) = uout_efcn(end);
    
    
else
    disp(['Iterations completed before continuation parameter reached. Current parameter value: ' num2str(par.(contPar.Name)) ' Desired value: ' num2str(contPar.final)])
end


figure(1); plot(cont_vec,free_vec1,'-o','LineWidth',2);
xlabel(contPar.Name); ylabel(contPar.Free1); 
title('Continuation Parameter');  box on;
set(gca,'fontsize',18,'linewidth',2);

figure(2); plot(cont_vec,free_vec2,'-o','LineWidth',2);
xlabel(contPar.Name); ylabel(contPar.Free2); 
title('Continuation Parameter'); box on;
set(gca,'fontsize',18,'linewidth',2);

plot_spiral(uout_sp,par,numPar);

% save final data
U = uout_sp(1:end-1);
Ubar = uout_efcn(1:end-1);


save(file_names.out_name ,'U','Ubar','par','numPar','cont_vec','free_vec1','free_vec2','cont_efcn');


