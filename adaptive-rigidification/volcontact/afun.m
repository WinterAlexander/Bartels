function [Ax] = afun( x, mesh, h )
% afun Evaluates ((1 - h * alpha0)M - (h alpha1 + h^2) K) x

mesh.dp = x;
df = computeForceDifferentials( mesh );
Ax = (1 - mesh.alpha0*h) * mesh.M * x - (h*mesh.alpha1 + h*h) * df;
