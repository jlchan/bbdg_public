#ifndef _FEM_HEADER
#define _FEM_HEADER

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dfloat.h"

#include "Mesh.h"
#include "Basis.h"

/* prototypes for storage functions (Utils.c) */
dfloat **BuildMatrix(int Nrows, int Ncols);
dfloat  *BuildVector(int Nrows);
int    **BuildIntMatrix(int Nrows, int Ncols);
int     *BuildIntVector(int Nrows);

dfloat **DestroyMatrix(dfloat **);
dfloat  *DestroyVector(dfloat *);
int    **DestroyIntMatrix(int **);
int     *DestroyIntVector(int *);

void PrintMatrix(char *message, dfloat **A, int Nrows, int Ncols);
void SaveMatrix(char *filename, dfloat **A, int Nrows, int Ncols);

/* geometric/mesh functions */
Mesh *ReadGmsh3d(char *filename);

void PrintMesh ( Mesh *mesh );

void Normals3d(Mesh *mesh, int k,
	       double *nx, double *ny, double *nz, double *sJ);


void GeometricFactors3d(Mesh *mesh, int k,
			double *drdx, double *dsdx, double *dtdx,
			double *drdy, double *dsdy, double *dtdy,
			double *drdz, double *dsdz, double *dtdz,
			double *J);

// start up
void StartUp3d(Mesh *mesh);
void BuildMaps3d(Mesh *mesh);
void FacePair3d(Mesh *mesh);

void InitQuadratureArrays(Mesh *mesh);

// set initial condition
void WaveSetU0(Mesh *mesh, dfloat *Q, dfloat time,
	       double(*uexptr)(double,double,double,double));
void WaveProjectU0(Mesh *mesh, dfloat *Q, dfloat time,
		   double(*uexptr)(double,double,double,double));

void WaveSetData3d(dfloat *Q);
void WaveGetData3d(Mesh *mesh, dfloat *Q);

// solver
void test_kernels(Mesh *mesh);
void time_kernels(Mesh *mesh);
void Wave_RK(Mesh *mesh, dfloat FinalTime, dfloat dt);

// for BB vs NDG stability (paper result)
void Wave_RK_sample_error(Mesh *mesh, dfloat FinalTime, dfloat dt,
                          double(*uexptr)(double,double,double,double));

void RK_step(Mesh *mesh, dfloat rka, dfloat rkb, dfloat fdt);
void compute_error(Mesh *mesh, double time, dfloat *Q,
		   double(*uexptr)(double,double,double,double),
		   double &L2err, double &relL2err);

void setOccaArray(MatrixXd A, occa::memory &B); // assumes matrix is double
void setOccaIntArray(MatrixXi A, occa::memory &B); // assumes matrix is int

void writeVisToGMSH(string fileName,Mesh *mesh, dfloat *Q, int iField, int Nfields);

// currently unused
void BuildFaceNodeMaps(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf,
		       MatrixXi &mapP);

#endif
