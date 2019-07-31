function [f,J]  = jacobian_for_spectra_karma_non_reflecting_BC(u,L1,L2,par,numPar)
% Solve Karma system with natural boundary conditions
% bndry_idx is the very outer radius
% Kappa is treated as a fixed parameter


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

tauE = 1/par.tauE;
taun = 1/par.taun;
R = 1/(1 - exp(-par.Re));
kappa = par.kappa;
omega = par.omega;

% Variables: Put into short grid format - only one point at the origin.
U = u(1:nx*ny); U = U(nx:end);
V = u(nx*ny+1:2*nx*ny); V = V(nx:end);

% Set up the boundary conditions
% Take radial derivative and select just the boundary 
Ur = Lr*U; %Ur = Ur(bndry_idx);
Vr = Lr*V; %Vr = Vr(bndry_idx);

% Select boundary, then take angular derivative
Uphi = Ltheta* U(bndry_idx);
Vphi = Ltheta* V(bndry_idx);

% Reaction terms
[thSn,thSn_prime] = thetaS(U - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(U - par.Eh);        
fE = tauE .* ( -U + 0.5*(par.Estar - V.^par.M).*(1 - tanE) .* U.^2 );
fn = taun .* ( R.* thSn - V);       

line1 = par.gamma.*(L2*U) + omega.*(L1*U) + fE;
line2 = par.delta.*(L2*V) + omega.*(L1*V) + fn;

% Boundary conditions
line1(bndry_idx) = 0; line2(bndry_idx) = 0;
line1(bndry_idx) = Ur - kappa.*Uphi;
line2(bndry_idx) = Vr - kappa.*Vphi;

f = [line1; line2];% Not really used, just for completeness
   
   
% Jacobian              
if (nargout > 1)
   
    u_bndry_idx = nx*(ny-2)+2: nx*(ny-1) + 1;
    v_bndry_idx = nx*(ny-1)+nx*(ny-2)+2 :  2*nx*(ny-1) + 1;
    
    % Nonlinear terms
    fE_E = tauE*(-1 - 0.5.*(par.Estar - V.^par.M).*(sech(par.Eh - U).^2).* (U.^2) + (par.Estar - V.^par.M).*(1 - tanE).*U);
    fE_n = tauE*(-0.5* par.M.*V.^(par.M-1) .* (1 - tanE).*(U.^2));
    fn_E = taun*(R.*thSn_prime);
    fn_n = -taun;
     
    I = speye(nx*(ny-1)+1,nx*(ny-1)+1);
    Ltheta_jac = zeros(size(Lr));
    Ltheta_jac(:,nx*(ny-2)+2:nx*(ny-1)+1) = Ltheta;
    
    dU = [par.gamma.*L2 + omega.*L1 + spdiags(fE_E ,0,nx*(ny-1)+1,nx*(ny-1)+1);
            spdiags(fn_E ,0,nx*(ny-1)+1,nx*(ny-1)+1)];
 
    dU(u_bndry_idx,:) = 0; dU(v_bndry_idx,:) = 0; 
    dU(u_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
     
    dV = [spdiags(fE_n,0,nx*(ny-1)+1,nx*(ny-1)+1);
       par.delta*L2 + omega*L1 + fn_n*I];
    
    dV(u_bndry_idx,:) = 0; dV(v_bndry_idx,:) = 0; 
    dV(v_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
       
    J = [dU, dV];
    J = sparse(J);
   
          
end


    
    
    
    
    
    
    