function [f] = computeForces( mesh )
% computeForces comptues element forces
%   computeForces( mesh, mu, lambda ) sets the force accumulators in the
%   mesh given the current state and given deformation parameters.

mu = mesh.mu;
lambda = mesh.lambda;
els = mesh.el;
p = mesh.p;
f = mesh.f;
for m = 1:size(els,1)
    el = els(m);
    i = el.t(1);
    j = el.t(2);
    k = el.t(3);
    pi = p(i*2-1:i*2);
    pj = p(j*2-1:j*2);
    pk = p(k*2-1:k*2);    
    Ds = [ pi-pk, pj-pk ];
    F = Ds * el.Bm;
    E = 1/2 * ( F' * F - eye(2) );
    P = F * ( 2 * mu * E + lambda * trace(E) * eye(2) );
    H = - el.area * P * el.Bm';
    f(i*2-1:i*2) = f(i*2-1:i*2) + H(:,1);
    f(j*2-1:j*2) = f(j*2-1:j*2) + H(:,2);
    f(k*2-1:k*2) = f(k*2-1:k*2) - H(:,1) - H(:,2);
end