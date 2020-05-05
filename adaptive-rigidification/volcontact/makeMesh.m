function [ mesh ] = makeMesh( p, t, rho,  mu, lambda, alpha0, alpha1 )
% makeMesh creates a FEM mesh

% TODO: might want to repack everything into vectors rather than nx2
% matrices
N = size(p,1);
mesh.N = N;                 % number of nodes
mesh.p = reshape(p',N*2,1); % position state
mesh.p0 = mesh.p;           % initial positions
mesh.v = zeros(N*2,1);      % velocity state
mesh.v0 = zeros(N*2,1);     % initial velocities
mesh.f = zeros(N*2,1);      % force accumulator
mesh.dp = zeros(N*2,1);     % for computing mesh differentials
mesh.df = zeros(N*2,1);     % force differential accumulator
mesh.t = t;                 % triangles
mesh.el = makeElements( p, t );
mesh.edges = makeBoundaries( p, t );
mesh.pinned = zeros(N,1);   % flags pinned indices
mesh.pinnedInds = [];       % list of pinned node IDs
mesh.pinnedDOFs = [];
mesh.unpinnedDOFs = 1:2*N;
mesh.mass = computeMass( mesh, rho ); % density of rubber
mesh.M = sparse( 1:2*N, 1:2*N, mesh.mass ); % sparse matrix form for convenience
mesh.mu = mu;           % Lamé parameter
mesh.lambda = lambda;   % Lamé parameter
mesh.alpha0 = alpha0;   % Rayleigh parameter
mesh.alpha1 = alpha1;   % Rayleigh parameter
mesh.DOFIndexOffset = 0; % offset of these DOF indexes in total list of indexes