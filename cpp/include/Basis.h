#ifndef _BASIS_HEADER
#define _BASIS_HEADER

#include "dfloat.h"

// testing eigen
void test_solve();
void test_basis();

// 1D 
void JacobiGQ(int N, int alpha_int, int beta_int, VectorXd &r, VectorXd &w);
VectorXd JacobiP(VectorXd x, dfloat alpha, dfloat beta, int d);
VectorXd GradJacobiP(VectorXd x, double alpha, double beta, int p);
MatrixXd Bern1D(int N, VectorXd r);

// 2D
void rstoab(VectorXd r, VectorXd s, VectorXd &a,VectorXd &b);
VectorXd Simplex2DP(VectorXd a,VectorXd b,int i, int j);
MatrixXd Vandermonde2D(int N, VectorXd r, VectorXd s); // for face
MatrixXd BernTri(int N, VectorXd r, VectorXd s);

// 3D
void rsttoabc(VectorXd r,  VectorXd s,  VectorXd t, 
	      VectorXd &a, VectorXd &b, VectorXd &c);
VectorXd Simplex3DP(VectorXd a,VectorXd b,VectorXd c,int i, int j, int k);
void GradSimplex3DP(VectorXd a, VectorXd b, VectorXd c,
		    int id, int jd, int kd,
		    VectorXd &V3Dr, VectorXd &V3Ds, VectorXd &V3Dt);
MatrixXd Vandermonde3D(int N, VectorXd r, VectorXd s, VectorXd t);
void GradVandermonde3D(int N, VectorXd r, VectorXd s, VectorXd t,
		       MatrixXd &V3Dr, MatrixXd &V3Ds, MatrixXd &V3Dt);

// unused routine - can define WB nodes this way
MatrixXd VandermondeGHsurf(int N, int Npsurf, VectorXd r, VectorXd s, VectorXd t);
// Gordon-Hall blending VDM for a single face
void VandermondeHier(int N, VectorXd r, VectorXd s, VectorXd t,
			 MatrixXd &vertV, MatrixXd &edgeV, MatrixXd &faceV);
void barytors(VectorXd L1, VectorXd L2, VectorXd L3, VectorXd &r,VectorXd &s);

// Bernstein-Bezier basis
MatrixXd BernTet(int N, VectorXd r, VectorXd s, VectorXd t);		 
void GradBernTet(int N,  VectorXd r, VectorXd s, VectorXd t,
		 MatrixXd &V1, MatrixXd &V2, MatrixXd &V3, MatrixXd &V4);

// geofacs
MatrixXd vgeofacs3d(VectorXd x,VectorXd y,VectorXd z,
		    MatrixXd Dr,MatrixXd Ds,MatrixXd Dt);
MatrixXd sgeofacs3d(VectorXd x,VectorXd y,VectorXd z,
		    MatrixXd Drf,MatrixXd Dsf,MatrixXd Dtf);

// extract tabulated nodes
void Nodes3D(int N, VectorXd &r, VectorXd &s, VectorXd &t);
void tet_cubature(int N, VectorXd &rq, VectorXd &sq, VectorXd &tq, VectorXd &wq);
void tri_cubature(int N, VectorXd &rq, VectorXd &sq, VectorXd &wq);
//void tet_cubature_duffy(int N, VectorXd &a, VectorXd &wa,
//			VectorXd &b,  VectorXd &wb,
//			VectorXd &c,  VectorXd &wc);

// helper routines not present in C++ std lib
unsigned int factorial(int n);
unsigned int factorial_ratio(int n1,int n2);
unsigned int nchoosek(int n,int k);

// linalg helper routines (emulate Matlab)
MatrixXd mldivide(MatrixXd &A, MatrixXd &B); // A\B
MatrixXd mrdivide(MatrixXd &A, MatrixXd &B); // A/B
VectorXd extract(const VectorXd& full, const VectorXi &ind); // currently unused
void get_sparse_ids(MatrixXd A, MatrixXi &cols, MatrixXd &vals); // get fixed-bandwidth sparse ids

#endif
