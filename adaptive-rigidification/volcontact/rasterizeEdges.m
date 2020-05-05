function [] = rasterizeEdges( mesh )
% rasterizeEdges Rasterize the edges of a mesh.

range = [-3, 3];
width = range(2)-range(1);
N = 256;
layers = zeros(N,maxLayers);

for i = 1:size(mesh.edges,1)
    p1x = mesh.p( mesh.edges(i,1), 1 );
    p1y = mesh.p( mesh.edges(i,1), 2 );
    p2x = mesh.p( mesh.edges(i,2), 1 );
    p2y = mesh.p( mesh.edges(i,2), 2 );

    % fit an affine function
    % evaluate on a range
    
    % TODO finish this thing!
    
end
