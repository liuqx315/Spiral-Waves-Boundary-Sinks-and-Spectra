function [F,J] = Karma_boundary_sink_jacobian(w,par,numPar,mesh_params)
% Jacobian for Karma model

global U0;

% Rename parameters
L1y = mesh_params.DT;
L2x = mesh_params.D2Z;
nt = mesh_params.nt;
nz = mesh_params.nz; nz_long = mesh_params.nz_long;
L2X_long = mesh_params.L2X_long;

wU = w(1:nt*nz);
wV = w(nt*nz+1:2*nz*nt);

tauE = 1/par.tauE;
taun = 1/par.taun;
gamma = par.gamma;
delta = par.delta;
R = 1/(1 - exp(-par.Re));


par.omega = 2*pi/par.Tperiod;
[u_ff0,U0,~,~] = get_karma_rolls(U0,par,numPar,mesh_params,mesh_params.Dt,mesh_params.Dt2);  
u_ff = u_ff0(1:nt*nz_long); v_ff = u_ff0(nt*nz_long+1:2*nt*nz_long);


dt = -2*pi./par.Tperiod; % Temporal scaling
dL = (par.kappa./(2*pi*mesh_params.N))^2;

% Reaction terms
fU = @(u,v) tauE .* ( -u + 0.5*(par.Estar - v.^par.M).*(1 - tanh(u-par.Eh)) .* u.^2 );
fV = @(u,v) taun .* ( R.* 0.5.*(1 + tanh(par.s.*(u - par.En))) - v);

D2u = dL*L2X_long * (mesh_params.chi_ffLong.*u_ff); D2u = D2u(mesh_params.bc_idx);
D2v = dL*L2X_long * (mesh_params.chi_ffLong.*v_ff); D2v = D2v(mesh_params.bc_idx);

u_ff_short = mesh_params.chi_ff.*u_ff(mesh_params.bc_idx);
v_ff_short = mesh_params.chi_ff.*v_ff(mesh_params.bc_idx);

line1 = dt.*L1y*(u_ff_short + wU) + gamma.*D2u + gamma.* (dL*L2x * wU) + fU(u_ff_short + wU, v_ff_short + wV); 
line2 = dt.*L1y*(v_ff_short + wV) + delta.*D2v + delta.* (dL*L2x * wV) + fV(u_ff_short + wU, v_ff_short + wV);

line1(mesh_params.iend) = wU(mesh_params.iend); % dirchlet bcs at the LHS of the domain for w
line2(mesh_params.iend) = wV(mesh_params.iend);

F = [line1; line2]; 

% Jacobain
if nargout > 1
        
    fE_E = @(u,v) tauE*(-1 - 0.5.*(par.Estar - v.^par.M).*(sech(par.Eh - u).^2).* (u.^2) + (par.Estar - v.^par.M).*(1 - tanh(u - par.Eh)).*u);
    fE_n = @(u,v) tauE*(-0.5* par.M.*v.^(par.M-1) .* (1 - tanh(u - par.Eh)).*(u.^2));
    fn_E = @(u,v) taun*(R.*0.5 .* par.s .* (sech(par.s.*(u - par.En)).^2));
    fn_n = -taun;
    
    
    I = speye(nt*nz,nt*nz);
    
   
     dwU = [dt.*L1y + gamma.*dL*L2x + spdiags(fE_E(u_ff_short + wU,v_ff_short+wV),0,nz*nt,nz*nt);
         spdiags(fn_E(u_ff_short+wU, v_ff_short + wV),0,nz*nt,nz*nt)];

     dwV = [spdiags(fE_n(u_ff_short + wU,v_ff_short+wV),0,nz*nt,nz*nt);
         dt.*L1y + delta.*dL*L2x + fn_n.*I];   
     
   
     J = [dwU, dwV];
   
     % Boundary conditions in Jacobian
     J(mesh_params.iend,:) = 0; J(nz*nt+mesh_params.iend, :) = 0; 
     J(mesh_params.iend,mesh_params.iend)=speye(length(mesh_params.iend));
     J(nz*nt+mesh_params.iend, nz*nt+mesh_params.iend)=speye(length(mesh_params.iend));    
     
     J = sparse(J);
     
end







