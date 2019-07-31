function [F,J] = Rossler_boundary_sink_jacobian(w,par,numPar,mesh_params)
% Solve the Rossler model on a half-infinite channel 

global U0;

% Rename parameters
L1y = mesh_params.DT;
L2x = mesh_params.D2Z;
nt = mesh_params.nt;
nz = mesh_params.nz;nz_long = mesh_params.nz_long;
L2X_long = mesh_params.L2X_long;

a = par.a;
b = par.b;
c = par.c;
delta1 = par.delta1;
delta2 = par.delta2;
delta3 = par.delta3;

wU = w(1:nt*nz);
wV = w(nt*nz+1:2*nz*nt);
wW = w(2*nt*nz+1:3*nz*nt);


[u_ff0,U0,~,~] = get_rossler_rolls(U0,par,numPar,mesh_params,mesh_params.Dt,mesh_params.Dt2);  
u_ff = u_ff0(1:nt*nz_long); v_ff = u_ff0(nt*nz_long+1:2*nt*nz_long);
w_ff = u_ff0(2*nt*nz_long+1:3*nt*nz_long);


dt = -2*pi./par.Tperiod; % Temporal scaling
dL = (par.kappa./(2*pi*mesh_params.N))^2; % Spatial scaling

% Reaction terms
fU = @(u,v,w) -v - w;
fV = @(u,v,w) u + a.*v;
fW = @(u,v,w) u.*w - c.*w + b;

D2u = dL*L2X_long * (mesh_params.chi_ffLong.*u_ff); D2u = D2u(mesh_params.bc_idx);
D2v = dL*L2X_long * (mesh_params.chi_ffLong.*v_ff); D2v = D2v(mesh_params.bc_idx);
D2w = dL*L2X_long * (mesh_params.chi_ffLong.*w_ff); D2w = D2w(mesh_params.bc_idx);

u_ff = mesh_params.chi_ff.*u_ff(mesh_params.bc_idx);  % Put on short grid
v_ff = mesh_params.chi_ff.*v_ff(mesh_params.bc_idx);
w_ff = mesh_params.chi_ff.*w_ff(mesh_params.bc_idx);
 
line1 = dt.*L1y*(u_ff + wU) + delta1.*D2u + delta1.* (dL*L2x * wU) + fU(u_ff + wU, v_ff + wV, w_ff + wW); 
line2 = dt.*L1y*(v_ff + wV) + delta2.*D2v + delta2.* (dL*L2x * wV) + fV(u_ff + wU, v_ff + wV, w_ff + wW);
line3 = dt.*L1y*(w_ff + wW) + delta3.*D2w + delta3.* (dL*L2x * wW) + fW(u_ff + wU, v_ff + wV, w_ff + wW);

line1(mesh_params.iend) = wU(mesh_params.iend); % dirchlet bcs at the LHS of the domain for w
line2(mesh_params.iend) = wV(mesh_params.iend);
line3(mesh_params.iend) = wW(mesh_params.iend);

F = [line1; line2; line3];

% Jacobain
if nargout > 1
    
    fU_V = -1;  fU_W = -1;
    fV_U = 1;   fV_V = a;
    fW_U = @(u,v,w) w;   fW_W = @(u,v,w) u - c;
    
    I = speye(nt*nz,nt*nz);
    
   
     dwU = [dt.*L1y + delta1.*dL*L2x ;
         fV_U.*I;
         spdiags(fW_U(u_ff + wU, v_ff + wV, w_ff + wW),0,nz*nt,nz*nt)];
    
     
     dwV = [fU_V.*I;
         dt.*L1y + delta2.*dL*L2x + fV_V.*I;
         sparse(nz*nt,nz*nt)];
     
     dwW = [fU_W.*I;
         sparse(nt*nz,nt*nz);
         dt.*L1y + delta3.*dL*L2x + spdiags(fW_W(u_ff + wU, v_ff + wV, w_ff + wW),0,nz*nt,nz*nt)];
    
     
     J = [dwU, dwV, dwW];
   
     % Boundary conditions in Jacobian
%      J(mesh_params.iend,:) = 0; J(nz*nt+mesh_params.iend, :) = 0; J(2*nz*nt+mesh_params.iend, :) = 0;
     J(mesh_params.iend,mesh_params.iend)=speye(length(mesh_params.iend));
     J(nz*nt+mesh_params.iend, nz*nt+mesh_params.iend)=speye(length(mesh_params.iend));
     J(2*nz*nt+mesh_params.iend, 2*nz*nt+mesh_params.iend)=speye(length(mesh_params.iend));     
%      
     J = sparse(J);
     
end







