function [m] = computeMass( mesh, rho )
% computeMass comptues the diagonal of the mass matrix
%   computeMass( mesh, rho ) returns the lumped node diagonal of the mass 
%   matrix for the given mesh

m = zeros( mesh.N*2,1 );
el = mesh.el;
for i = 1:size(el,1)
    for j = 1:3
        m( el(i).t(j)*2-1 ) = m( el(i).t(j)*2-1 ) + 1/3 * el(i).area * rho;
        m( el(i).t(j)*2 ) = m( el(i).t(j)*2 ) + 1/3 * el(i).area * rho;
    end
end