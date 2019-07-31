%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerically continue absolute spectrum of spiral wave 
% Assumes you have an initial starting point for the absolute spectrum
% Uses the 2pi-periodic wave train of the spiral wave
% Continues along spectral curve and computes all spatial eigenvalues to
% ensure overlapping spatial eigenvalue criteria for absolute spectrum is met
% Simple continuation
% Stephanie Dodson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

%% Select system to solve
% Modify for Karma or Rossler
file_names.starting_point = 'absolute_spectrum_data/karma_abs_spec_re0p6_start.mat';  % Starting point on the absolute spectrum
file_names.out_name = 'absolute_spectrum_data/karma_abs_spec_re0p6_branch1_plus.mat';

file_names.cont_func = 'absolute_spectra_karma';   % Function that continues along the absolute spectrum
file_names.spatial_evals_fcn = 'spatial_evals_fcn_karma';  % Function that computes the spatial eigenvalues (independent of continuation)


%%
soln = load(file_names.starting_point);
u_infty = soln.u_infty;                  % Asymptotic wave trains
numPar = soln.numPar;
par = soln.par;


U0 = soln.uout(1:par.numVars*numPar.nx);              % Eigenvector corresponding to par.nu
W0 = zeros(size(U0));              % Second eigenvectors (of second spatial eigenvalue par.nu + 1i*par.beta) - random guesses, which will be solved for


contPar.ds = 0.1;                         % Continuation step size: sign determines direction of continuation
contPar.Name = 'beta';                      % Continuation parameter: vertical separation between spatial eigenvalues: Im(nu_1) - Im(nu_1)
contPar.free1 = 'nu';                   % Free parameter 1
contPar.free2 = 'lambda';                       % Free parameter 2
contPar.numContSteps =100;
contPar.final = 0;                       % Stopping condition
contPar.sp_evals_iter = 5;              % Computes spatial eigenvalues and updates plots ever contPar.sp_evals_iter iterations

% Additional variables
options = optimset('Display','iter','Jacobian','on', 'DerivativeCheck','off',...
                    'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',5000);

addpath ../Util/
fh = str2func(file_names.cont_func);
fh_spatial = str2func(file_names.spatial_evals_fcn);
%% Initial data

% Compute linear operator
[L1, L2] = ComputeLinearOperator_1D(par,numPar); % 1st and 2nd derivative matrices: assumes 2pi-periodic domain

% Initial condition is the loaded data
initial_sol = [U0; W0; par.(contPar.free1); par.(contPar.free2)]; % Eigenvectors
phase_cond.u_old = U0;
phase_cond.w_old = W0;

% Do an fsolve to make sure eigenvector is right
uout = fsolve(@(y) fh(y,u_infty,L1,L2,par,numPar,phase_cond),initial_sol,options);
par.(contPar.free1) = uout(end-1);
par.(contPar.free2) = uout(end);

phase_cond.u_old = uout(1:par.numVars*numPar.nx);
phase_cond.w_old = uout(par.numVars*numPar.nx+1:2*par.numVars*numPar.nx);

%% Continuation

free1_vec = zeros(contPar.numContSteps+1,1);   % temporal eigenvalues
cont_vec = zeros(contPar.numContSteps+1,1);       % gamma
free2_vec = zeros(contPar.numContSteps+1,1);       % Spatial eigenvalues
diff_vals = zeros(contPar.numContSteps,1);

free1_vec(1) = par.(contPar.free1);
free2_vec(1) = par.(contPar.free2);
cont_vec(1) = par.(contPar.Name);

for k = 1:contPar.numContSteps
    disp(k)
    
    par.(contPar.Name) = par.(contPar.Name) + contPar.ds; % Update continuation parameter: simple continuation
    
    disp([contPar.Name ': ' num2str(par.(contPar.Name))])
    
    uout = fsolve(@(y) fh(y,u_infty,L1,L2,par,numPar,phase_cond),uout,options);
    par.(contPar.free1) = uout(end-1);
    par.(contPar.free2) = uout(end);
    phase_cond.u_old = uout(1:par.numVars*numPar.nx);
    phase_cond.w_old = uout(par.numVars*numPar.nx+1:2*par.numVars*numPar.nx);
    
    free1_vec(k+1) = par.(contPar.free1);
    free2_vec(k+1) = par.(contPar.free2);
    cont_vec(k+1) = par.(contPar.Name);

    % Find spatial eigenvalues - mostly for debugging
    
    if mod(k,contPar.sp_evals_iter) == 0
        
        [A,vals] = fh_spatial(L1,L2,u_infty,par,numPar);

            figure(1); hold('on'); plot(vals,'.','MarkerSize',28); title(['Spatial Eigenvalues, \lambda = ' num2str(par.lambda) ', \nu = ' num2str(par.nu)]);
            plot(par.nu,'s','MarkerSize',15,'MarkerEdgeColor','black','MarkerFaceColor','black');   % black square denotes par.nu - spatial eigenvalue used in the continuation
            xlabel('Re(\nu)'); ylabel('Im(\nu)'); box on; set(gca,'fontsize',14,'linewidth',2);  axis([-1,1,-10,10]); drawnow;  % axis will need to be adjusted for different models
             
            figure(2); plot(free2_vec(1:k),'-o'); title('Temporal Eigenvalues','fontsize',16); 
            xlabel('Re(\lambda)'); ylabel('Im(\lambda)'); box on; set(gca,'FontSize',18,'linewidth',2); drawnow;

    end
    
        % Additional stopping criteria
       if par.(contPar.Name)*sign(contPar.ds) >= contPar.final*sign(contPar.ds)
           disp('Continuation Parameter Reached');
           
           free1_vec = free1_vec(1:k+1);
           free2_vec = free2_vec(1:k+1);
          
           break
       end
    

end



%save(file_names.out_name ,'uout','par','numPar','free1_vec','cont_vec','u_infty','free2_vec','contPar');

