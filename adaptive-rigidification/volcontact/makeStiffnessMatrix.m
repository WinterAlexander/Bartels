function [ K ] = makeStiffnessMatrix( mesh )
% makeStiffnessMatrix builds a sparse stiffness matrix
%
%    Building the stiffness matrix is generally a bad idea, and this method
%    for doing so is not particularly fast.

N = mesh.N;
K = sparse( N*2, N*2 );
for i = 1:2*N
    mesh.dp = zeros( N*2, 1 );
    mesh.dp(i) = 1;
    K(:,i) = computeForceDifferentials( mesh );
end