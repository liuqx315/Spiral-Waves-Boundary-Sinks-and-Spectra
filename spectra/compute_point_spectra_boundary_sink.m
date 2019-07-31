% Calculate the eigenvalues of the boundary sink
close all; clear all;


%% Select the system
% Modify for Rossler or Karma
file_names.boundary_sink_file = '../data_files/Karma_boundary_sink_re0p6.mat';
file_names.out_file = '../data_files/Karma_boundary_sink_re0p6_pointspectra.mat';

file_names.contruct_initial_solution = 'get_karma_rolls';   % Construct the far-field solution and initial condition
file_names.jacobian = 'Karma_boundary_sink_jacobian';

%% Set up
load(file_names.boundary_sink_file);

% Set the seed points for eigS computation - computed numVals near each seed point
numVals = 10;

% Define the real and imaginary values for seed points
seed_real = -10:5:0;       
seed_imag = [-abs(par.omega),0,abs(par.omega)];

[LR,LI] = meshgrid(seed_real,seed_imag);  % Make a grid of seeds
LR = LR(:); LI = LI(:);
loc = LR + 1i.*LI;  % Grid of seed points

addpath ../Util/

fh_rolls = str2func(file_names.contruct_initial_solution);
fh_jacobian =  str2func(file_names.jacobian);

plot_full_boundary_sink_pattern(wU,U0,par,numPar,mesh_params,fh_rolls); % Plot the full boundary sink pattern: far-field + boundary region.

%% Calculate the eigenvalues

% Find the jacobain
[F,J] = fh_jacobian([wU;par.Tperiod], par,numPar,mesh_params);

vals = cell(length(loc),1);
vecs = cell(length(loc),1);

disp(['Total: ' num2str(length(loc))]);

figure; hold on;
for j = 1:length(loc)
    disp(j)
    [vecs{j},tmp] = eigs(J,numVals,loc(j));
    vals{j} = diag(tmp);
    
    plot(real(vals{j}),imag(vals{j}),'k.','MarkerSize',18); drawnow;
    
end

plot([0,0],ylim,'k-','linewidth',2)
plot(xlim,[0,0],'k-','linewidth',2)
xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
set(gca,'fontsize',20);

save(file_names.out_file,'vals','vecs','par','numPar','U0', 'wU','mesh_params','-v7.3')
