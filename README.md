## BBDG: a Bernstein-Bezier Discontinuous Galerkin solver for CPUs or GPUs
## Copyright (c) 2016, Jesse Chan, Tim Warburton
## https://github.com/jlchan/bbdg_public
## BBDG is released under the GPL lesser public license. see license.bbdg for details.

BBDG is based on the Matlab codes of "Nodal Discontinuous Galerkin Methods" by Jan Hesthaven and Tim Warburton, as well as the MIDG code of Tim Warburton. The implementation is based on the approach outlined in arxiv.org/abs/1512.06025. 

BBDG is built with Eigen 3.2.8 (most recent stable version as of Feb 16, 2016), which is included in this repository. 
To update Eigen, replace header files in include/Eigen with new header files. (http://eigen.tuxfamily.org)

BBDG depends on the OCCA library (http://libocca.org) for portability between different GPU vendors. 

High order visualization can be achieved using GMSH (http://gmsh.info).

These codes have been tested on OSX 10.11.16.

- To compile with a degree 3 nodal basis 

> `make N=3 B=0 -j`

- To compile with a degree 3 Bernstein-Bezier basis:

> `make N=3 B=1 -j`

- To run:

> `./main meshes/cube1.msh`


