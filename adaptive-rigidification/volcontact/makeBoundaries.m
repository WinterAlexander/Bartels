function [edges] = makeBoundaries( p, t )
% makeBoundaries Creates a list of boundareis from triangles.
%   makeBoundaries( p, t ) returns an array of boundary edges given 2D
%   points p (stored as rows), and triangles t defined by indices
%   (1 indexed and stored as rows).
% 
%   Assumes there is less than 100000 points in the mesh (gross yes).

N = 100000;
if ( size(p,1) >= N ) 
    error('triangles must have less than 100000 points');
end
s = java.util.HashSet();
for i = 1:size(t,1)
    for j = 1:3
        a = t(i,j);
        b = t(i,mod(j,3)+1);
        hash1 = a*N + b;
        hash2 = b*N + a;
        if ( s.contains( hash2 ) ) 
            s.remove( hash2 );
        else
            s.add( hash1 );
        end
    end
end
L=s.toArray;
edges = zeros(size(L,1),2);
for i = 1:size(L,1)
    edges(i,:) = [floor(L(i)/N), mod(L(i),N)];
end
