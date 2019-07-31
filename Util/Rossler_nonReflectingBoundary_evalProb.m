function [f,J]  = Rossler_nonReflectingBoundary_evalProb(u,u_sp, par,numPar, phase_cond)
% Solve Rossler system with natural boundary conditions
% Solve eigenvalue problem as well
% bndry_idx is the very outer radius
% Kappa is treated as a fixed parameter

L1 = numPar.L1;
L2 = numPar.L2./(par.r2^2);
Lr = numPar.Lr./(par.r2);
Ltheta = numPar.Ltheta;
bndry_idx = numPar.bndry_idx;

nx = numPar.nx;
ny = numPar.ny;

a = par.a;
b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;
kappa = par.kappa;
omega = par.omega;

U = u_sp(1:nx*ny);
V = u_sp(nx*ny+1:2*nx*ny);
W = u_sp(2*nx*ny+1:3*nx*ny);

Ubar = u(1:nx*ny);
Vbar = u(nx*ny+1:2*nx*ny);
Wbar = u(2*nx*ny+1:3*nx*ny); 
lambda = u(end);


% Set up the boundary conditions
% Take radial derivative and select just the boundary 
UrBar = Lr*Ubar;
VrBar = Lr*Vbar;
WrBar = Lr*Wbar;

% Select boundary, then take angular derivative
UphiBar = Ltheta* Ubar(bndry_idx);
VphiBar = Ltheta* Vbar(bndry_idx);
WphiBar = Ltheta* Wbar(bndry_idx);

% Reaction terms
% fU = -V - W;
% fV = U + a.*V;
% fW = U.*W - c.*W + b;

% Linearizations
fU_V = -1;  fU_W = -1;
fV_U = 1;   fV_V = a;
fW_U = W;   fW_W = U - c;

% Solve eigenvalue problem
line1 = delta1.*(L2*Ubar) + omega.*(L1*Ubar) + fU_V.*Vbar + fU_W.*Wbar - lambda.*Ubar;
line2 = delta2.*(L2*Vbar) + omega.*(L1*Vbar) + fV_U.*Ubar + fV_V.*Vbar - lambda.*Vbar;
line3 = delta3.*(L2*Wbar) + omega.*(L1*Wbar) + fW_U.*Ubar + fW_W.*Wbar - lambda.*Wbar;
line4 = phase_cond.Ubarold'*Ubar - 1;  % Eigenfunction normalization

% Boundary conditions
line1(bndry_idx) = 0; line2(bndry_idx) = 0; line3(bndry_idx) = 0;
line1(bndry_idx) = UrBar - kappa.*UphiBar;
line2(bndry_idx) = VrBar - kappa.*VphiBar;
line3(bndry_idx) = WrBar - kappa.*WphiBar;

f = [line1; line2; line3; line4];
   
% Jacobian              
if (nargout > 1)
   
    u_bndry_idx = nx*(ny-1)+1:nx*ny; v_bndry_idx = nx*ny+nx*(ny-1)+1:2*nx*ny; w_bndry_idx = 2*nx*ny+nx*(ny-1)+1:3*nx*ny;
    I = speye(nx*ny,nx*ny); Z = sparse(nx*ny,nx*ny);
    Ltheta_jac = zeros(size(Lr)); Ltheta_jac(:,nx*(ny-1)+1:nx*ny) = Ltheta;
        
    dUbar = [delta1.*L2 + omega.*L1 - lambda.*I;
        fV_U.*I; spdiags(fW_U,0,nx*ny,nx*ny);  phase_cond.Ubarold' ];
    dUbar(u_bndry_idx,:) = 0; dUbar(v_bndry_idx,:) = 0; dUbar(w_bndry_idx,:) = 0;
    dUbar(u_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    dVbar = [ fU_V.*I; delta2.*L2 + omega.*L1 + (fV_V - lambda).*I; 
        Z; sparse(1,nx*ny) ]; 
    dVbar(u_bndry_idx,:) = 0; dVbar(v_bndry_idx,:) = 0; dVbar(w_bndry_idx,:) = 0;
    dVbar(v_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    dWbar = [fU_W.*I; Z; delta3.*L2 + omega.*L1 + spdiags(fW_W - lambda,0,nx*ny,nx*ny) ; sparse(1,nx*ny) ]; 
    dWbar(u_bndry_idx,:) = 0; dWbar(v_bndry_idx,:) = 0; dWbar(w_bndry_idx,:) = 0;
    dWbar(w_bndry_idx,:) = Lr - kappa.*Ltheta_jac;

    dlambda = [ -Ubar; -Vbar; -Wbar; 0];
    dlambda(u_bndry_idx) = 0; dlambda(v_bndry_idx) = 0; dlambda(w_bndry_idx) = 0;
  
    J = [ dUbar, dVbar, dWbar, dlambda];
    J = sparse(J);
   
          
end


    
    
    
    
    
