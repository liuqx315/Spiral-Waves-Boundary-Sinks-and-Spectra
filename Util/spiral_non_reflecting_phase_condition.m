% Set up differentiation matrices and phase condition for the
% non-reflecting boundary conditions

% Differentiation matrices
[L1, L2] = ComputeLinearOperator(par,numPar);       % 2D differentiation matrices: d/dtheta and radial laplacian
[Ltheta,~] = ComputeLinearOperator_1D(par,numPar);  % 1D theta differentiation matrix

% Radial derivative
h = 1/(numPar.ny-1);
Dr = zeros(1,numPar.ny);
Dr(end-4:end) = [1/4; -4/3; 3; -4; 25/12]./h;

Ix = speye(numPar.nx,numPar.nx);
Lr = kron(Dr,Ix);

numPar.bndry_idx = (numPar.nx*(numPar.ny-1)+1):numPar.nx*numPar.ny;  % indices for boundary 

% Store matrices in structure
numPar.L1 = L1; numPar.L2 = L2;
numPar.Ltheta = Ltheta; numPar.Lr = Lr;


% Phase condition
pc = zeros(numPar.nx,numPar.ny);
pc(:,ceil(numPar.ny/2)) = 1; 
phase_cond.pc = pc(:);%reshape(pc,numPar.nx*numPar.ny,1);

us = (phase_cond.pc).*U0(1:numPar.nx*numPar.ny);
phase_cond.u_star = us(phase_cond.pc>0);

us_th = (L1*U0(1:numPar.nx*numPar.ny)); % d/d(theta) of phase condition
phase_cond.u_star_th = us_th(phase_cond.pc > 0);
