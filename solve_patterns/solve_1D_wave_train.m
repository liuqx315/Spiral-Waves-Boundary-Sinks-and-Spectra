%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically solve for a 2pi-periodic wave train in moving coordinate
% frame
% Assumes you have an initial guess
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; 

%% System to solve
% Modify for Karma or Rossler
file_names.wave_train = '../data_files/Rossler_1D_wave_train_c2.mat';         % Initial guess
file_names.out_file = '../data_files/Rossler_1D_wave_train_c2_solved.mat';    % Name of output saved data

file_names.wt_solver = 'Rossler_1D';      % Function to solve the 1D wave train

%% Set up
u0 = load(file_names.wave_train);
U0 = u0.U;
numPar = u0.numPar;
par = u0.par;


% Set the free parameter
par.Free = 'omega';

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',... % fsolve options
                   'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',1500);
      
addpath ../Util/

fh = str2func(file_names.wt_solver);
%% Do an fsolve to find solution
[L1, L2] = ComputeLinearOperator_1D(par,numPar);
initial_sol = [U0; par.(par.Free)];
phase_cond.u_old = U0(1:numPar.nx);

uout = fsolve(@(y) fh(y,L1,L2,par,numPar,phase_cond),initial_sol,options);  
U = uout(1:par.numVars*numPar.nx);
par.(par.Free) = uout(end);

save(file_names.out_file,'par','numPar','U');



