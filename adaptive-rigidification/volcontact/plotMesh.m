function [] = plotMesh( mesh )
% plotMesh plots given edges, mesh, and pinned vertices (assumes hold on)

pr = reshape( mesh.p, 2, mesh.N )';
edges = mesh.edges;

patch('vertices',pr,'faces',mesh.t,'edgecol',[.5,.5,.5],'facecol',[.75,.75,.75]);
    
plot( [pr(edges(:,1),1)'; pr(edges(:,2),1)'], [pr(edges(:,1),2)'; pr(edges(:,2),2)'], 'k');

plot( pr(mesh.pinnedInds,1), pr(mesh.pinnedInds,2), 'r.' );
