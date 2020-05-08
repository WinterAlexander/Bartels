% create some triangle meshes with http://persson.berkeley.edu/distmesh/

%addpath('../distmesh')


draw = 0;
skip = 1;
N = 450;        % number of steps
g = -9.8;       % gravity

rho = 100;%1522;     % density of rubber
nu = 0.48;      % Poisson ratio: close to 0.5 for rubber
k = 2e4;     % Young's modulus: 0.01e9 approximate for rubber
mu = k / 2 / (1 - nu);                  % Lamé parameter
lambda = k * nu / (1+nu) / (1-2*nu);    % Lamé parameter
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K

h = 0.01;%0.01; %0.00003;     % time step

SE = 0; % Symplectic Euler or Backward Euler
tol = 1e-6;  
maxit = 15;

clear meshes;

% Good for a single triangle test:
%p = [ 0, 0; 1 0; 1 1 ];
%t = [ 1 2 3 ];
%meshes(1) = makeMesh( p, t, rho, mu, lambda );
%meshes(1).pinnedInds = [1,2];% find(p(:,2)<-2.49 );
%meshes(1).pinned( meshes(1).pinnedInds ) = 1;
%%meshes(1).p(3*2-1) = meshes(1).p(3*2-1) + 0.2;
%meshes(1).v(3*2-1) = meshes(1).v(3*2-1) + 100;

fd=@(p) drectangle(p,-1.5,1.5,-2.5,-1.5);
[p,t]=distmesh2d(fd,@huniform,0.3,[-1.5,-2.5;1.5,-1.5],[-1.5,-2.5;1.5,-2.5;-1.5,-1.5;1.5,-1.5]);
meshes(1) = makeMesh( p, t, rho, mu, lambda, alpha0, alpha1 );
meshes(1).pinnedInds = find(p(:,2)<-2.49 );
%meshes(1).pinnedInds = find(p(:,1)<-1.49 );
meshes(1).pinned( meshes(1).pinnedInds ) = 1;
%meshes(1).v(3*2-1) = meshes(1).v(3*2-1) + 0.4;  % initial perturbation
meshes(1).pinnedDOFs = sort( [ meshes(1).pinnedInds *2 - 1; meshes(1).pinnedInds*2 ] );
meshes(1).unpinnedDOFs = setdiff( 1:size(meshes(1).f,1), meshes(1).pinnedDOFs );
meshes(1).DOFIndexOffset = 0;

fd=@(p) sqrt(sum(p.^2,2))-1;
[p,t]=distmesh2d(fd,@huniform,0.3,[-1,-1;1,1],[]);
meshes(2) = makeMesh( p, t, rho, mu, lambda, alpha0, alpha1 );
meshes(2).DOFIndexOffset = meshes(1).DOFIndexOffset + meshes(1).N*2;

elapsed = 0;


if draw == 1
    v = VideoWriter('anim4.avi');
    open(v);
end

for i = 1:size(meshes, 2)
    V = reshape(meshes(i).p, size(meshes(i).p, 1)/2, 2);
    meshes(i).B = linear_tri2dmesh_B(V, meshes(i).t);
end

for i = 1:N

    % Draw the meshes
    if (draw == 1 && mod(i,skip) == 0 )
        %subplot(2,1,1);  % uncomment to make other plots in the same figure for debugging!
        cla;
        hold on;
        for j = 1:size(meshes,2)
           plotMesh( meshes(j) );
        end    
        axis equal;
        axis([-3,3,-4,2]);
        title(sprintf('t = %f', elapsed ));
    end

    %subplot(2,1,2);
    [V, dVdp] = rayTraceEdges( meshes );
    axis equal;
    axis([-3,3,-4,2]);
    
    if draw == 1
        F = getframe();
        writeVideo(v,F);
    end
    
    
    % Advance the state of the meshes
    for j = 1:size(meshes,2)    
        mesh = meshes(j);

        if ( SE ) 
            mesh.f = zeros( mesh.N*2, 1 );
            % Gravity
            mesh.f(2:2:end) = mesh.mass(2:2:end) * g;
            % Elastic forces
            mesh.f = computeForces( mesh );
            % Rayleigh damping
            mesh.f = mesh.f - mesh.alpha0 * mesh.M * mesh.v;        
            mesh.dp = mesh.v;
            mesh.df = computeForceDifferentials( mesh );
            mesh.f = mesh.f + mesh.alpha1 * mesh.df; % note signs
            % pressure force
            mesh.f = mesh.f - dVdp( mesh.DOFIndexOffset + (1:mesh.N*2))';

            mesh.v = mesh.v + h * diag(1./mesh.mass) * mesh.f; 
            mesh.p = mesh.p + h * mesh.v;   
            mesh.p(mesh.pinnedDOFs) = mesh.p0(mesh.pinnedDOFs); 
        else
            % compute right hand side
            rhs = zeros( mesh.N*2, 1 );
            % Gravity
            % Could split the ODE if poorly conditioned to make sure that 
            % stuff still falls under gravity 
            %hg = zeros( mesh.N*2, 1 ); 
            %hg(2:2:end) = h*g;
            rhs(2:2:end) = h* mesh.mass(2:2:end) * g;
            % Elastic forces
            rhs = rhs + h *computeForces( mesh );
            % Rayleigh damping and stifness linearization (implicitization)               
            mesh.dp = mesh.v;
            mesh.df = computeForceDifferentials( mesh );
            rhs = rhs + (h^2 + h*mesh.alpha1) * mesh.df; % note signs
            rhs = rhs + h * mesh.alpha0 * mesh.M * mesh.v;
            rhs = rhs - h * 10000 * dVdp( mesh.DOFIndexOffset + (1:mesh.N*2))';

            
            % Need to remove DOFs...  :(  This would be way faster for an
            % already slow matlab implementation
            %fcn = @(x) afun(x, mesh, h);
            %[deltav, fl, rr, it, rv] = pcg( fcn, rhs, tol, maxit );
            
            % SLOTH!!!
          
            
            F = mesh.B * mesh.p;
            
            for n = 1:size(F, 1)
                if mod(n, 9) == 0
                    F(n) = 1;
                end
            end
            
            params = zeros(size(mesh.t, 1), 2);
            for n = 1:size(mesh.t, 1)
                params(n, 1:2) = [mesh.mu, mesh.lambda];
            end
            
            C = d2psi_neohookean_dF2(mesh.t, F, params);
            K = -B * C * B; 
            A = ((1-h*mesh.alpha0)*mesh.M-(h*mesh.alpha1+h^2)*K);
            
            ii = mesh.unpinnedDOFs;
            deltav = zeros( 2*mesh.N, 1 );
            deltav(ii) = A(ii,ii) \ rhs(ii);
            
%            if (j == 1) 
%                subplot(2,1,2);
%                plot( rhs );
%                title('debug deltav vector');
%            end
                
            mesh.v = mesh.v + deltav;% + hg'; 
            mesh.p = mesh.p + h * mesh.v;   
            mesh.p(mesh.pinnedDOFs) = mesh.p0(mesh.pinnedDOFs); 
        end
        
        meshes(j) = mesh;
    end
    elapsed = elapsed + h;
end

close(v);

% TODO: 
% - Experiment with different mesh resolutions (and shapes?)
% - Try SE=1, adjust the steping and stiffness to not explode
% - Add the second mesh
% - Fix backward euler

% - compute contacts with LDI approach
% - start with volume based penalty
% - extend to contact constraints with friction
% - extend from mono-volume to multi-volume
