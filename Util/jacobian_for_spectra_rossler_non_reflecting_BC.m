function [f,J] = jacobian_for_spectra_rossler_non_reflecting_BC(u,L1,L2,par,numPar)
% J = Jacobian on the reduced polar grid, with only one grid point at the
% origin
% Natural boundary conditions
% Need  L1, L2 in short grid format

nx = numPar.nx;
ny = numPar.ny;

% Extra matrices needed for non-reflecting boundary conditions
[Ltheta,~] = ComputeLinearOperator_1D(par,numPar);  % 1D theta differentiation matrix


h = par.r2/(ny-1);
Dr = zeros(1,ny-1);
Dr(end-4:end) = [1/4; -4/3; 3; -4; 25/12]./h;

Ix = speye(nx,nx);
Lr = kron(Dr,Ix);
Lr = [zeros(nx,1), Lr];
Lr = sparse(Lr);

bndry_idx = (nx*(ny-2)+2) : ( nx*(ny-1) + 1);

a = par.a;
b = par.b;
c = par.c;
omega = par.omega;
kappa = par.kappa;

delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);
W = u(2*nx*ny+1:3*nx*ny); W = W(nx:end);

% Set up the boundary conditions
% Take radial derivative and select just the boundary 
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

line1(bndry_idx) = 0; line2(bndry_idx) = 0; line3(bndry_idx) = 0;
line1(bndry_idx) = Ur - kappa.*Uphi;
line2(bndry_idx) = Vr - kappa.*Vphi;
line3(bndry_idx) = Wr - kappa.*Wphi;

f = [line1; line2; line3]; % Not really used, just for completeness
   
% Jacobian              
if (nargout > 1)
    
    u_bndry_idx = nx*(ny-2)+2: nx*(ny-1) + 1;
    v_bndry_idx = nx*(ny-1)+nx*(ny-2)+2 :  2*nx*(ny-1) + 1;
    w_bndry_idx = 2*nx*(ny-1)+nx*(ny-2)+2 :  3*nx*(ny-1) + 1;
    
    fU_V = -1;  fU_W = -1;
    fV_U = 1;   fV_V = a;
    fW_U = W;   fW_W = U - c;
    
    I = speye(nx*(ny-1)+1,nx*(ny-1)+1);
    Ltheta_jac = zeros(size(Lr));
    Ltheta_jac(:,nx*(ny-2)+2:nx*(ny-1)+1) = Ltheta;
        
    dU = [delta1.*L2 + omega.*L1;
        fV_U.*I;
        spdiags(fW_U,0,nx*(ny-1)+1,nx*(ny-1)+1)];
    dU(u_bndry_idx,:) = 0; dU(v_bndry_idx,:) = 0; dU(w_bndry_idx,:) = 0;
    dU(u_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    dV = [fU_V.*I;
        delta2.*L2 + omega.*L1 + fV_V.*I;
        sparse(nx*(ny-1)+1,nx*(ny-1)+1)];
    dV(u_bndry_idx,:) = 0; dV(v_bndry_idx,:) = 0; dV(w_bndry_idx,:) = 0;
    dV(v_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
    
    dW = [fU_W.*I;
        sparse(nx*(ny-1)+1,nx*(ny-1)+1);
        delta3.*L2 + omega.*L1 + spdiags(fW_W,0,nx*(ny-1)+1,nx*(ny-1)+1)];
    
    dW(u_bndry_idx,:) = 0; dW(v_bndry_idx,:) = 0; dW(w_bndry_idx,:) = 0;
    dW(w_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
    
        
    J = [dU, dV, dW];
    J = sparse(J);
   
          
end




