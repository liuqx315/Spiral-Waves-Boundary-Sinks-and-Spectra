function [f,J]  = Rossler_2D_spiral_nonreflecting(u,par,numPar, phase_cond)
% Solve Rossler system with natural boundary conditions
% bndry_idx is the very outer radius
% Kappa is treated as a fixed parameter
% Free parameter: omega

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

U = u(1:nx*ny);
V = u(nx*ny+1:2*nx*ny);
W = u(2*nx*ny+1:3*nx*ny);
omega = u(end); 

% Select the solution evaluated at R/2
u_phase = U.*phase_cond.pc;
u_phase = u_phase(phase_cond.pc>0); 

% Set up the boundary conditions
% Take radial derivative, which selects just the boundary 
Ur = Lr*U; 
Vr = Lr*V; 
Wr = Lr*W; 

% Select boundary, then take angular derivative
Uphi = Ltheta* U(bndry_idx);
Vphi = Ltheta* V(bndry_idx);
Wphi = Ltheta* W(bndry_idx);

% Reaction terms
fU = -V - W;
fV = U + a*V;
fW = U.*W - c.*W + b;

line1 = delta1.*(L2*U) + omega.*(L1*U) + fU;
line2 = delta2.*(L2*V) + omega.*(L1*V) + fV;
line3 = delta3.*(L2*W) + omega.*(L1*W) + fW;
line4 = (phase_cond.u_star_th)'*(u_phase - phase_cond.u_star);

% Boundary conditions
line1(bndry_idx) = 0; line2(bndry_idx) = 0; line3(bndry_idx) = 0;
line1(bndry_idx) = Ur - kappa.*Uphi;
line2(bndry_idx) = Vr - kappa.*Vphi;
line3(bndry_idx) = Wr - kappa.*Wphi;

f = [line1; line2; line3; line4];
   
% Jacobian              
if (nargout > 1)
   
    u_bndry_idx = nx*(ny-1)+1:nx*ny;
    v_bndry_idx = nx*ny+nx*(ny-1)+1:2*nx*ny;
    w_bndry_idx = 2*nx*ny+nx*(ny-1)+1:3*nx*ny;
    
    fU_V = -1;  fU_W = -1;
    fV_U = 1;   fV_V = a;
    fW_U = W;   fW_W = U - c;
    
    I = speye(nx*ny,nx*ny);
    phase_jacob = [sparse(1,numPar.nx*(ceil(numPar.ny/2) - 1)), phase_cond.u_star_th', sparse(1,numPar.nx*floor(numPar.ny/2)) ]; % Only the middle radii
    
    Ltheta_jac = zeros(size(Lr));
    Ltheta_jac(:,nx*(ny-1)+1:nx*ny) = Ltheta;
    
    
    dU = [delta1.*L2 + omega.*L1;
        fV_U.*I;
        spdiags(fW_U,0,nx*ny,nx*ny);
        phase_jacob];
 
    dU(u_bndry_idx,:) = 0; dU(v_bndry_idx,:) = 0; dU(w_bndry_idx,:) = 0;
    dU(u_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    
    dV = [fU_V.*I;
        delta2.*L2 + omega.*L1 + fV_V.*I;
        sparse(nx*ny,nx*ny);
        sparse(1,nx*ny)];
    
    dV(u_bndry_idx,:) = 0; dV(v_bndry_idx,:) = 0; dV(w_bndry_idx,:) = 0;
    dV(v_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    dW = [fU_W.*I;
        sparse(nx*ny,nx*ny);
        delta3.*L2 + omega.*L1 + spdiags(fW_W,0,nx*ny,nx*ny);
        sparse(1,nx*ny)];
    
    dW(u_bndry_idx,:) = 0; dW(v_bndry_idx,:) = 0; dW(w_bndry_idx,:) = 0;
    dW(w_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    domega1 = L1*U; domega1(bndry_idx) = 0;
    domega2 = L1*V; domega2(bndry_idx) = 0;
    domega3 = L1*W; domega3(bndry_idx) = 0;
    
    domega = [domega1; domega2; domega3; 0];
  
    J = [dU, dV, dW, domega];
    J = sparse(J);
   
          
end


    
    
    
    
    
