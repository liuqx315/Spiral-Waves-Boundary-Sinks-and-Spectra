%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find point eigenvalues of pre-computed spiral wave
% Uses sparse methods (eigs) with seed points defined by 'seed_real' and
% 'seed_imag'
% Works for neumann or non-reflecting boundary conditions, depending on
% Jacobian function
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;

%% Select system to solve
% Modify for Karma or Rossler and non-reflecting or neumann boundary
% conditions
file_names.spiral_file = '../data_files/Karma_spiral_R5_re0p6.mat';   % Precomputed spiral file
file_names.out_file = 'point_spectrum/Karma_spiral_R5_re0p6_ptSpec1.mat';

file_names.jacobian_fcn = 'jacobian_for_spectra_karma';  % Function to compute Jacobian: uses short grid. Function defines Neumann or Non-reflecting boundary conditions and Rossler or Karma

%% Set up
load(file_names.spiral_file);

% Set the seed points for eigS computation - computed numVals near each seed point
numVals = 10;

% Define the real and imaginary values for seed points
seed_real = -10:5:0;       
seed_imag = [-abs(par.omega),0,abs(par.omega)];

[LR,LI] = meshgrid(seed_real,seed_imag);  % Make a grid of seeds
LR = LR(:); LI = LI(:);
loc = LR + 1i.*LI;  % Grid of seed points

addpath ../Util/

[L1, L2] = ComputeLinearOperator_shortGrid(par,numPar); % This function does not use domain size of 1

fh_jacobian = str2func(file_names.jacobian_fcn);

[f,J] = fh_jacobian(U,L1,L2,par,numPar);


%% Spectra

vals = cell(length(loc),1);
vecs = cell(length(loc),1);

disp(['Total: ' num2str(length(loc))]);

figure; hold on;

for j = 1:length(loc)
    disp(j)
    [vecs{j},tmp] = eigs(J,numVals,loc(j));
    vals{j} = diag(tmp);
    
    plot(real(vals{j}),imag(vals{j}),'.','MarkerSize',18); drawnow;
    
end


xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
box on; set(gca,'fontsize',20,'linewidth',2);

save(out_file,'vals','par','numPar','vecs')



