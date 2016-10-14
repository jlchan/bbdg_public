#ifndef __DFLOAT_DEFINED

#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
#include <stdio.h>

#define dfloat float
#define dfloat4 float4
//#define dfloat double
//#define dfloat4 double4

//#define MPI_DFLOAT MPI_FLOAT

// dfloat matrix with eigen
#define MatrixXdf Matrix<dfloat,Dynamic,Dynamic>

#define __DFLOAT_DEFINED
#endif
