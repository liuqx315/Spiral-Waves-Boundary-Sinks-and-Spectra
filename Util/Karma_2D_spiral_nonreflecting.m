function [f,J]  = Karma_2D_spiral_nonreflecting(u,par,numPar, phase_cond)
% Solve Karma system with non-reflecting boundary conditions
% bndry_idx is the very outer radius
% Kappa is treated as a fixed parameter
% Free parameter: omega

% Rename parameters and differentiation matrices
L1 = numPar.L1;
L2 = numPar.L2./(par.r2^2);
Lr = numPar.Lr./(par.r2);
Ltheta = numPar.Ltheta;
bndry_idx = numPar.bndry_idx;

nx = numPar.nx;
ny = numPar.ny;

tauE = 1/par.tauE;
taun = 1/par.taun;
R = 1/(1 - exp(-par.Re));
kappa = par.kappa;

U = u(1:nx*ny);
V = u(nx*ny+1:2*nx*ny);
omega = u(end); 

% Select the solution evaluated at R/2
u_phase = U.*phase_cond.pc;
u_phase = u_phase(phase_cond.pc>0); 

% Set up the boundary conditions
% Take radial derivative, which selects just the boundary 
Ur = Lr*U; 
Vr = Lr*V; 

% Select boundary, then take angular derivative
Uphi = Ltheta* U(bndry_idx);
Vphi = Ltheta* V(bndry_idx);

% Reaction terms
[thSn,thSn_prime] = thetaS(U - par.En,par.s);   % Smoothed Heaviside function
tanE = tanh(U - par.Eh);        % Only calculate once
fE = tauE .* ( -U + 0.5*(par.Estar - V.^par.M).*(1 - tanE) .* U.^2 );
fn = taun .* ( R.* thSn - V);       % No origin smoothing

line1 = par.gamma.*(L2*U) + omega.*(L1*U) + fE;
line2 = par.delta.*(L2*V) + omega.*(L1*V) + fn;
line3 = (phase_cond.u_star_th)'*(u_phase - phase_cond.u_star);

% Boundary conditions
line1(bndry_idx) = 0; line2(bndry_idx) = 0;
line1(bndry_idx) = Ur - kappa.*Uphi;
line2(bndry_idx) = Vr - kappa.*Vphi;

f = [line1; line2; line3];
   
% Jacobian              
if (nargout > 1)
   
    u_bndry_idx = nx*(ny-1)+1:nx*ny;
    v_bndry_idx = nx*ny+nx*(ny-1)+1:2*nx*ny;
    
    % Nonlinear terms
    fE_E = tauE*(-1 - 0.5.*(par.Estar - V.^par.M).*(sech(par.Eh - U).^2).* (U.^2) + (par.Estar - V.^par.M).*(1 - tanE).*U);
    fE_n = tauE*(-0.5* par.M.*V.^(par.M-1) .* (1 - tanE).*(U.^2));
    fn_E = taun*(R.*thSn_prime);
    fn_n = -taun;
    
    phase_jacob = [sparse(1,numPar.nx*(ceil(numPar.ny/2) - 1)), phase_cond.u_star_th', sparse(1,numPar.nx*floor(numPar.ny/2)) ]; % Only the middle radii
    
    Ltheta_jac = zeros(size(Lr));
    Ltheta_jac(:,nx*(ny-1)+1:nx*ny) = Ltheta;
    
    dU = [par.gamma.*L2 + omega.*L1 + spdiags(fE_E ,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny);
            spdiags(fn_E ,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny) ;
            phase_jacob];
 
    dU(u_bndry_idx,:) = 0; dU(v_bndry_idx,:) = 0; 
    dU(u_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
     
    dV = [spdiags(fE_n,0,numPar.nx*numPar.ny,numPar.nx*numPar.ny);
       par.delta*L2 + omega*L1 + fn_n*speye(numPar.nx*numPar.ny);  % No smoothing function
       sparse(1,numPar.nx*numPar.ny)];
    
    dV(u_bndry_idx,:) = 0; dV(v_bndry_idx,:) = 0; 
    dV(v_bndry_idx,:) = Lr - kappa.*Ltheta_jac;
        
    domega1 = L1*U; domega1(bndry_idx) = 0;
    domega2 = L1*V; domega2(bndry_idx) = 0;
    domega = [domega1; domega2; 0];
    
    J = [dU, dV, domega];
    J = sparse(J);
   
          
end


    
    
    
    
    
    
    