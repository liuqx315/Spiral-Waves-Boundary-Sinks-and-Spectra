function [u_ff0,U,Du_ff,wU0] = get_rossler_rolls(u,par,numPar,mesh_params,L1,L2)
% Solves the 1D wave train problem with
% Input: u - [U;V] periodic wave train
%        L1, L2: 1D periodic differentiation matrices (for wave train
%        problem)
% Output: u_ff0, v_ff0: wave train expanded to the (x,t)-plane long grid
%         U: solved wave train
% .       Du_ff: derivative of wave train


phase_cond.u_old = u(1:numPar.nx);
u1 = [u; par.omega];


options = optimset('Display','off','Jacobian','on', 'DerivativeCheck','off',...
                    'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',200);
                
               
uout = fsolve(@(y) Rossler_1D(y,L1,L2,par,numPar,phase_cond),u1,options);
U = uout(1:numPar.nx);
V = uout(numPar.nx+1:2*numPar.nx);
W = uout(2*numPar.nx+1:3*numPar.nx);


%% Interpolate onto x-t mesh

u_ff0 = zeros(mesh_params.nt,mesh_params.nz_long);
v_ff0 = zeros(mesh_params.nt,mesh_params.nz_long);
w_ff0 = zeros(mesh_params.nt,mesh_params.nz_long);
Du_ff = zeros(mesh_params.nt,mesh_params.nz_long);

Du = L1*U;

tmpU = repmat(U,mesh_params.N + mesh_params.Nc,1)';
tmpV = repmat(V,mesh_params.N + mesh_params.Nc,1)';
tmpW = repmat(W,mesh_params.N + mesh_params.Nc,1)';
tmpDU = repmat(Du,mesh_params.N + mesh_params.Nc,1)';

for j = 1:numPar.nx
    
    u_ff0(j,:) = tmpU;
    v_ff0(j,:) = tmpV;
    w_ff0(j,:) = tmpW;
    Du_ff(j,:) = tmpDU;
    
     tmpU = [tmpU(2:end), tmpU(1)];  % Shift based on the direction of group velocity
    tmpV = [tmpV(2:end), tmpV(1)]; 
    tmpW = [tmpW(2:end), tmpW(1)];
    tmpDU = [tmpDU(2:end), tmpDU(1)];
    

end

% Put into vector form for next use
U = [U;V;W];
u_ff0 = u_ff0(:); v_ff0 = v_ff0(:); w_ff0 = w_ff0(:);
Du_ff = Du_ff(:);

wU0 = mesh_params.chi_p.*u_ff0(mesh_params.bc_idx); % initial condition for the localized function 
wV0 = mesh_params.chi_p.*v_ff0(mesh_params.bc_idx);
wW0 = mesh_params.chi_p.*w_ff0(mesh_params.bc_idx);

u_ff0 = [u_ff0; v_ff0; w_ff0];

wU0 = [wU0; wV0; wW0];

