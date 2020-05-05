function [df] = computeForceDifferentials( mesh )
% computeForceDifferentials comptues element force differentials

mu = mesh.mu;
lambda = mesh.lambda;
els = mesh.el;
p = mesh.p;
dp = mesh.dp;
df = zeros(mesh.N*2,1);
for m = 1:size(els,1)
    el = els(m);
    i = el.t(1);
    j = el.t(2);
    k = el.t(3);
    pi = p(i*2-1:i*2);
    pj = p(j*2-1:j*2);
    pk = p(k*2-1:k*2);
    dpi = dp(i*2-1:i*2);
    dpj = dp(j*2-1:j*2);
    dpk = dp(k*2-1:k*2);
    Ds = [ pi-pk, pj-pk ];
    dDs = [ dpi-dpk, dpj-dpk ];    
    F = Ds * el.Bm;
    dF = dDs * el.Bm;
    E = 1/2 * ( F' * F - eye(2) );
    dE = 1/2 * ( dF' * F + F' * dF );
    %P = F * ( 2 * mu * E + lambda * trace(E) * eye(2) );
    dP = dF * ( 2 * mu * E + lambda * trace(E) * eye(2) ) + F * ( 2 * mu * dE + lambda * trace(dE) * eye(2) );
    %H = - el.area * P * el.Bm';
    dH = - el.area * dP * el.Bm';     
    df(i*2-1:i*2) = df(i*2-1:i*2) + dH(:,1);
    df(j*2-1:j*2) = df(j*2-1:j*2) + dH(:,2);
    df(k*2-1:k*2) = df(k*2-1:k*2) - dH(:,1) - dH(:,2);
end