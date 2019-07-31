function [F,J] = Rossler_boundary_sink(w,par,numPar,mesh_params)
% Solve the Rossler model on a half-infinite channel 

global U0;

% Rename parameters
L1y = mesh_params.DT;
L2x = mesh_params.D2Z;
nt = mesh_params.nt;
nz = mesh_params.nz; nz_long = mesh_params.nz_long;
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
Tperiod = w(end);


par.omega = -2*pi/Tperiod;
[u_ff0,U0,Du_ff,~] = get_rossler_rolls(U0,par,numPar,mesh_params,mesh_params.Dt,mesh_params.Dt2);  
u_ff = u_ff0(1:nt*nz_long); v_ff = u_ff0(nt*nz_long+1:2*nt*nz_long);
w_ff = u_ff0(2*nt*nz_long+1:3*nt*nz_long);

dt = -2*pi./Tperiod; % Temporal scaling
dL = (par.kappa./(2*pi*mesh_params.N))^2;

% Reaction terms
fU = @(u,v,w) -v - w;
fV = @(u,v,w) u + a.*v;
fW = @(u,v,w) u.*w - c.*w + b;

D2u = dL*L2X_long * (mesh_params.chi_ffLong.*u_ff); D2u = D2u(mesh_params.bc_idx);
D2v = dL*L2X_long * (mesh_params.chi_ffLong.*v_ff); D2v = D2v(mesh_params.bc_idx);
D2w = dL*L2X_long * (mesh_params.chi_ffLong.*w_ff); D2w = D2w(mesh_params.bc_idx);

u_ff_short = mesh_params.chi_ff.*u_ff(mesh_params.bc_idx);
v_ff_short = mesh_params.chi_ff.*v_ff(mesh_params.bc_idx);
w_ff_short = mesh_params.chi_ff.*w_ff(mesh_params.bc_idx);
Du_ff= mesh_params.chi_ff.*Du_ff(mesh_params.bc_idx);
 
line1 = dt.*L1y*(u_ff_short + wU) + delta1.*D2u + delta1.* (dL*L2x * wU) + fU(u_ff_short + wU, v_ff_short + wV, w_ff_short + wW); 
line2 = dt.*L1y*(v_ff_short + wV) + delta2.*D2v + delta2.* (dL*L2x * wV) + fV(u_ff_short + wU, v_ff_short + wV, w_ff_short + wW);
line3 = dt.*L1y*(w_ff_short + wW) + delta3.*D2w + delta3.* (dL*L2x * wW) + fW(u_ff_short + wU, v_ff_short + wV, w_ff_short + wW);

line1(mesh_params.iend) = wU(mesh_params.iend); % dirchlet bcs at the LHS of the domain for w
line2(mesh_params.iend) = wV(mesh_params.iend);
line3(mesh_params.iend) = wW(mesh_params.iend);

% Phase condition
L_cut = 1/mesh_params.N;
iBm = find( mesh_params.xx(:) <= L_cut ); % find indices close to the end of the first roll

z = linspace(0,mesh_params.Lz,mesh_params.nz);  hz= z(2) - z(1);
iim= find(z <= L_cut); nnzm = length(iim);

wz = [1, 2*ones(1,nnzm-2)+2*mod([1:nnzm-2],2),1]*hz/3;   % Simpson weights for intergration int = w*u
wt = 2*mesh_params.Lt*ones(mesh_params.nt,1)/mesh_params.nt;
wwm= kron(wz,wt');
wwm= wwm(:)';
  
wBm = wU(iBm);                       % find w on the domain x = 0:2*pi/kappa and y = 0..Ly
u_prime = Du_ff(iBm);
line4 = wwm*(u_prime .* wBm);

F = [line1; line2; line3; line4];  


% Jacobain
if nargout > 1
    
    fU_V = -1;  fU_W = -1;
    fV_U = 1;   fV_V = a;
    fW_U = @(u,v,w) w;   fW_W = @(u,v,w) u - c;
    
    I = speye(nt*nz,nt*nz);
    
   phase_jacob = zeros(1,nz*nt);
   phase_jacob(iBm) = wwm*spdiags(u_prime,0,nnzm*mesh_params.nt,nnzm*mesh_params.nt);
    
     dwU = [dt.*L1y + delta1.*dL*L2x ;
         fV_U.*I;
         spdiags(fW_U(u_ff_short + wU, v_ff_short + wV, w_ff_short + wW),0,nz*nt,nz*nt);
         phase_jacob];
    
     
     dwV = [fU_V.*I;
         dt.*L1y + delta2.*dL*L2x + fV_V.*I;
         sparse(nz*nt,nz*nt);
         sparse(1,nz*nt)];
     
     dwW = [fU_W.*I;
         sparse(nt*nz,nt*nz);
         dt.*L1y + delta3.*dL*L2x + spdiags(fW_W(u_ff_short + wU, v_ff_short + wV, w_ff_short + wW),0,nz*nt,nz*nt);
         sparse(1,nz*nt)];
     
     epsiF = 1e-8;
     dF = Rossler_boundary_sink([w(1:end-1); w(end)+epsiF],par,numPar,mesh_params);
     dT = (dF - F)./epsiF;  
     
     J = [dwU, dwV, dwW, dT];
   
     % Boundary conditions in Jacobian
     J(mesh_params.iend,:) = 0; J(nz*nt+mesh_params.iend, :) = 0; J(2*nz*nt+mesh_params.iend, :) = 0;
     J(mesh_params.iend,mesh_params.iend)=speye(length(mesh_params.iend));
     J(nz*nt+mesh_params.iend, nz*nt+mesh_params.iend)=speye(length(mesh_params.iend));
     J(2*nz*nt+mesh_params.iend, 2*nz*nt+mesh_params.iend)=speye(length(mesh_params.iend));     
     
     J = sparse(J);
     
end







