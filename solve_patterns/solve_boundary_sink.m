% Numerically solve for a boundary sink

close all; clear all;
global U0;

%% System to solve
%Modify for Karma/Rossler
file_names.wave_train = '../data_files/Karma_boundary_sink_re0p6.mat';  % Input file
file_names.out_file = 'Karma_boundary_sink_re0p6_solved.mat';    % Name of output file

file_names.contruct_initial_solution = 'get_karma_rolls';   % Construct the far-field solution and initial condition
file_names.boundary_problem = 'Karma_boundary_sink';        % Function to solve for the pattern

%% Set up
wt = load(file_names.wave_train);
par = wt.par; numPar = wt.numPar; % Corresponds to the wave train system and numerical parameters
u_wt = wt.U0;               % Wave train - used to compute initial condition for the boundary sink


% Additional parameters: Set up differentiation matrices and things
par.Tperiod = abs(2*pi/par.omega);  % Temporal period
par.Free = 'Tperiod';
N = 6;      % Number of rolls/wave trains in solution
Nc = 2;     % Number of wave trains to cut off LHS of final solution (remove boundary effects)

addpath ../Util/
set_up_matrices_boundary_sink;  

options = optimset('Jacobian','on','Display','iter','MaxIter',50,'DerivativeCheck','off');

fh_rolls = str2func(file_names.contruct_initial_solution);
fh_solver = str2func(file_names.boundary_problem);

%% Solve for pattern
% Form initial condition
% get far-field rolls
[~,U0,~,wU0] = fh_rolls(u_wt,par,numPar,mesh_params,Dt,Dt2);
plot_boundary_sink(wU0,mesh_params);

w0 = [wU0; par.Tperiod];

% Solve with fsolve
uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),w0,options);
wU = uout(1:end-1);
par.(par.Free) = uout(end);

% Plot results
plot_boundary_sink(wU,mesh_params);   % Plot part isolated at boundary
plot_full_boundary_sink_pattern(wU,u_wt,par,numPar,mesh_params,fh_rolls); % Plot the full boundary sink pattern: far-field + boundary region.

save(file_names.out_file,'wU','par','numPar','mesh_params','U0')




