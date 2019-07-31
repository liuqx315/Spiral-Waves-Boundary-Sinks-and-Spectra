%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically continue 1D 2pi periodic wave train in a parameter value
% Assumes you have an initial guess
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
%% Select system to solve
% Modify for Karma or Rossler
file_names.wave_train = '../data_files/Karma_1D_wave_train_re0p6.mat';  % contains initial starting solution/guess
file_names.out_file = '../data_files/Karma_1D_wave_train_re1p0.mat';    % Name of output saved data

file_names.wt_solver = 'Karma_1D';      % Function to solve the 1D wave train

%% Set up
% Load data
u0 = load(file_names.wave_train);
U0 = u0.U;              % Guess
numPar = u0.numPar;     % Structure of numerical parameters
par = u0.par;           % Structure of model parameters


% Continuation parameters
contPar.final = 1.0;
contPar.ds = 1;                % continuation step size: should always be positive!
contPar.initial_ds = 0.1;      % initial step: sign determines direction of continuation
contPar.Name = 'Re';            % Continuation parameter
contPar.Free = 'omega';         % Free parameter
contPar.numContSteps = 300;     % Maximum number of continuation steps
contPar.plot_iter = 20;         % Update bifurcation diagram/continuation plot ever plot_iter steps

% Additional variables
addpath ../Util/
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...  % fsolve options
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',500); 
fh = str2func(file_names.wt_solver);
x = numPar.Lx.*(0:numPar.nx-1)./numPar.nx;

%% Do an fsolve - ensure input data is a solution
[L1, L2] = ComputeLinearOperator_1D(par,numPar);  % 1st and 2nd derivative matrices: assumes 2pi-periodic domain!
initial_sol = [U0; par.(contPar.Free)];
phase_cond.u_old = U0(1:numPar.nx);

uout = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),initial_sol,options);  
par.(contPar.Free) = uout(end);

sol0 = [uout; par.(contPar.Name)];
phase_cond.u_old = uout(1:numPar.nx);

% Second point for secant method
par.(contPar.Name) = par.(contPar.Name) + contPar.initial_ds; % Update continuation parameter

uout = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),uout,options); 
par.(contPar.Free) = uout(end);
sol1 = [uout; par.(contPar.Name)]; % Second point
phase_cond.u_old = uout(1:numPar.nx);

cont_vec = zeros(1,contPar.numContSteps);
free_vec = zeros(1,contPar.numContSteps);
free_vec(1) = par.(contPar.Free);
cont_vec(1) = par.(contPar.Name);

for k = 1:contPar.numContSteps
    disp(k)
    disp([contPar.Name ': ' num2str(par.(contPar.Name))])
    
    % Predictor
    pred = sol1 + (sol1 - sol0)/norm(sol1 - sol0, 2) * contPar.ds;
    
    [uout,fval,exitflag,output,jacobian] = fsolve(@(sol) FixedPointSecantPredictorCorrector_1D_general(sol,sol1,sol0,par,numPar,L1,L2,contPar,phase_cond,fh),pred,options);
    par.(contPar.Free) = uout(end-1);
    par.(contPar.Name) = uout(end);
    
    phase_cond.u_old = uout(1:numPar.nx);

    cont_vec(k+1) = par.(contPar.Name);
    free_vec(k+1) = par.(contPar.Free);
   
    
    if mod(k,contPar.plot_iter) == 0
    
        figure(1); plot(cont_vec(1:k+1),free_vec(1:k+1),'-o'); title('Continuation Diagram','FontSize',16);
        ylabel(par.(contPar.Free),'FontSize',14); xlabel(contPar.Name,'FontSize',14); drawnow;
        
        
       
    end
    
    if (par.(contPar.Name)*sign(contPar.initial_ds)) >= (contPar.final*sign(contPar.initial_ds))

        disp(['Final Value reached: ' num2str(par.(contPar.Name))]);
        break;
    end
    
    % Prepare for next step
    sol0 = sol1; 
    sol1 = uout; 
    

end

if k < contPar.numContSteps
    % Set continuation parameter to exact value and do a final fsolve (if
    % final parameter reached)
    par.(contPar.Name) = contPar.final;
    uout = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),uout(1:end-1),options);  

    
else
    disp(['Continuation steps completed before final value reached. Final value: ' num2str(par.(contPar.Name))]);
    

end

U = uout(1:par.numVars*numPar.nx);
par.(contPar.Free) = uout(end);

save(file_names.out_file,'par','numPar','U','cont_vec','free_vec','contPar');



