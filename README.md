# meshMC
Mesh based Monte Carlo radiation transport


MCRT code that uses manifold triangular meshes(in .obj or .ply format) as the medium. Fluence is tracked on an underlying voxel structure.
Currently buggy.


## Monte Carlo on a Gourd
  ![MC Gourd](https://raw.githubusercontent.com/lewisfish/meshMC/master/gourdMC.png)

TODO:
- [ ] Implment kd-tree to accelerate photon triangle intersection
- [ ] Implment Fresnel reflections at mesh surfaces
- [ ] Ability to use multiple meshes
- [ ] Fix bugs
