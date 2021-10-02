# COMS

COMS is a conservative overset mesh solver. Overset mesh technique is useful in problems involving multiple relatively moving components. With overset mesh technique each component is meshed independently. Since there is no face-to-face connectivity between overset meshes, flux cannot flow through cell-faces across component meshes. This causes non-conservation of flux in mesh-system. COMS physically connects each component mesh by remeshing inter-grid regions with Advancing Front Delaunay algorithm. This way total flux in mesh-system is conserved.

COMS solves the two-dimensional Euler equations on unstructured grid topology. [Tailor](https://github.com/orxshi/tailor) which is a three-dimensional and load-balancing overset mesh solver is an improved version of COMS however, Tailor does not remesh component meshes into a single grid, therefore, is not flux-conservative.

Below is an example shows a background and an overset mesh. Some of the cells in the overlapping region are removed and remeshed to connect the meshes.

![](https://github.com/orxshi/COMS/blob/master/images/front.gif)

The following figure shows application of COMS to oscillating airfoil test case. At each time step, COMS remeshes the inter-grid region.

![](https://github.com/orxshi/COMS/blob/master/images/osc.gif)
