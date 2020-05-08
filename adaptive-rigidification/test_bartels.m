[V,T,F] = readMESH('data/meshes_mesh/coarser_bunny.mesh');

disp('Mesh loaded');

B = linear_tetmesh_dphi_dX(V, T)