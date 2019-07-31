
% Set up the partition functions, matrices, other numerical parameters for
% the boundary sink problem

% Spatial coordinates: x direction
nz = N*numPar.nx; Lz = 1; hz = Lz/(nz-1); z = linspace(0,Lz,nz);
nz_long = (N+Nc)*numPar.nx; 
tmp = (numPar.nx*hz).*(0:numPar.nx-1)./(numPar.nx-1); tmp = -flip(tmp);
z_long = linspace(-Nc*numPar.nx*hz,Lz,nz_long);
Lz_long = z_long(end) - z_long(1); 
hz_long = Lz_long/(nz_long-1); 

% Temporal discretisation in t
nt = numPar.nx; Lt = 2*pi; ht = 2*pi/nt;  t = ht*(1:nt); 


% Partition functions
d = 4*Lz/5; % cut off point at x = d
m = 40;  % stepness parameter
mesh_params.d = d;

% Differentiation matrix in z: finite differences 
ez = ones(nz,1); 
Dz = spdiags([ez -8*ez 0*ez 8*ez -ez],-2:2, nz, nz);
Dz(1,:)   = 0; Dz(2,2)   = 1;
Dz(nz,:) = 0; Dz(nz-1,nz-1) = -1;
Dz = Dz/(12*hz);        % first derivative, 4th order. Not sure on BC - think Neumann with ghost point

D2z = spdiags([-ez 16*ez -30.*ez 16*ez -ez], -2:2, nz, nz);
D2z(1,2)= 32; D2z(1,3)= -2;
D2z(2,1)= 16; D2z(2,2)= -31; 
D2z(nz,:) = D2z(1,nz:-1:1);
D2z(nz-1,:) = D2z(2,nz:-1:1);
D2z = D2z/(12*hz^2);  % second derivative, 4th order. Think Neumann BC with ghost point - has domain size in it!

% Longer grid
ez_long = ones(nz_long,1); 
L2X_long = spdiags([-ez_long 16*ez_long -30.*ez_long 16*ez_long -ez_long], -2:2, nz_long, nz_long);
L2X_long(1,2)= 32; L2X_long(1,3)= -2;
L2X_long(2,1)= 16; L2X_long(2,2)= -31; 
L2X_long(nz_long,:) = L2X_long(1,nz_long:-1:1);
L2X_long(nz_long-1,:) = L2X_long(2,nz_long:-1:1);
L2X_long = L2X_long/(12*hz_long^2);  % second derivative, 4th order. Think Neumann BC with ghost point - has domain size in it!


Iz = speye(nz);
wz = [1, 2*ones(1,nz-2)+2*mod([1:nz-2],2),1]/3;   % Simpson weights for intergration int = w*u un-scaled  must *hx

% Time derivatives
[~, Dt]  = fourdif(nt,1);  % 1st derivative matrix, Fourier
[~, Dt2] = fourdif(nt,2); % 2nd derivative matrix, Fourier: Need for the wave train

wt = 2*pi*ones(nt,1)/nt; % integration weights for trapzoid rule - mean
wzt = kron(wz,wt'); wzt=wzt(:)';

It = speye(nt); 

% Linear 2D differentiation matrices
DT    = kron(Iz,Dt);
DZ    = kron(Dz  ,It);
D2Z   = kron(D2z,It);
L2X_long = kron(L2X_long,It);

[xx,tt] = meshgrid(z,t); % 2D mesh
[xx2,tt2] = meshgrid(z_long,t); % Long grid mesh

% Cut-off functions
chi_p = 1/2 + 1/2*tanh(m*(xx(:)-d));  % 0 on left, 1 on right - not centered
chi_m = chi_p(end:-1:1);    % 1 on left, 0 on right - not centered
chi_ff = 1-chi_p;

chi_pLong = 1/2 + 1/2*tanh(m*(xx2(:)-d));  % 0 on left, 1 on right - not centered
chi_ffLong = 1-chi_pLong;

% Remove LHS boundary impacts on u_ff
% LHS Boundary conditions: want half-line, so will use Neumann, but on a
% larger domain than is actually solved

bc = ones(nt,nz_long);
bc(:,1:Nc*numPar.nx) = 0;
bc = bc(:);
bc_idx = find(bc>0);

iend = find(abs(xx)==0);

mesh_params.nz  = nz;     mesh_params.nt  = nt;
mesh_params.Lz  = Lz;     mesh_params.Lt  = Lt;
mesh_params.xx  = xx;     mesh_params.tt  = tt;
mesh_params.DT = DT;      mesh_params.D2Z = D2Z; mesh_params.DZ = DZ;
mesh_params.Dt = Dt;      mesh_params.Dt2 = Dt2;
mesh_params.xx2  = xx2;   mesh_params.tt2  = tt2;
mesh_params.N = N;        mesh_params.nz_long = nz_long;
mesh_params.Nc = Nc;      mesh_params.wzt = wzt;
mesh_params.L2X_long = L2X_long;
mesh_params.bc_idx = bc_idx;
mesh_params.iend=iend;
mesh_params.chi_ff = chi_ff;  
mesh_params.chi_ffLong = chi_ffLong;
mesh_params.chi_p = chi_p;


