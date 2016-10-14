#include "fem.h"

#include <iostream>
#include <fstream>

double *VX, *VY, *VZ;


Mesh *ReadGmsh3d(char *filename){
  
  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));

  std::ifstream meshfile(filename);
  if(!meshfile.is_open()){
    std::cout << "BBDG ERROR: Mesh '" << filename << "' not opened" << std::endl;
    exit(1);
  }

  int dummy, allVertices, allElements;

  while(meshfile.good()){

    std::string line;
    getline(meshfile, line);

    /// scan vertices in gmsh file
    if(line.find("$Nodes") != std::string::npos){

      /// read number of vertices
      meshfile >> allVertices;
      mesh->Nv = allVertices;

      /// resize vectors for the coordinates of vertices
      VX = (double*) calloc(mesh->Nv,sizeof(double));
      VY = (double*) calloc(mesh->Nv,sizeof(double));
      VZ = (double*) calloc(mesh->Nv,sizeof(double));

      /// scan line of each vertex, ignoring gmsh index of vertex
      for (int v = 0; v < allVertices; ++v){
        meshfile >> dummy >> VX[v] >> VY[v] >> VZ[v];
	//printf("VXYZ = [%.4f, %.4f, %.4f]\n",VX[v],VY[v],VZ[v]);
      }
    }

    /// scan elements in gmsh file
    if(line.find("$Elements") != std::string::npos){

      /// read number of elements
      meshfile >> allElements;
      getline(meshfile, line);

      /// record location in file
      int marker = meshfile.tellg();
      int eGmshNum, eGeoType, Tags;
      VectorXi EToGmshE(allElements);

      /// count number of triangular/tetrahedral elements
      int KTri = 0;
      int KTetra = 0;
      int extra=0;
      for (int k = 1;k <= allElements; ++k) {
        meshfile >> eGmshNum >> eGeoType;
        if(eGeoType == 2){  /// for triangular element
          KTri++;
	} else if(eGeoType == 4){  /// for tetrahedral element
	  EToGmshE(KTetra) = eGmshNum;
          KTetra++;
	}else{
          extra++;
	}
        getline(meshfile, line);

      }
      std::cout << extra << " unknown element types read in gmsh file" << std::endl;
      
      /// rewind to beginning of element list
      meshfile.seekg(marker);

      /// shape EToV to hold Elements triples
      mesh->Nverts = 4; /* assume tets */
      mesh->Nedges = 6; /* assume tets */
      mesh->Nfaces = 4; /* assume tets */

      mesh->K = KTetra;

      mesh->EToV = BuildIntMatrix(KTetra, mesh->Nverts);
      mesh->EToGmshE = EToGmshE;
      //mesh->GX = BuildMatrix(mesh->K, mesh->Nverts);
      //mesh->GY = BuildMatrix(mesh->K, mesh->Nverts);
      //mesh->GZ = BuildMatrix(mesh->K, mesh->Nverts);
      mesh->GX.resize(mesh->K, mesh->Nverts);
      mesh->GY.resize(mesh->K, mesh->Nverts);
      mesh->GZ.resize(mesh->K, mesh->Nverts);

      int *ETag = BuildIntVector(KTetra);
      int *ENumGmsh = BuildIntVector(KTetra);
      for (int k = 0; k < KTetra; ++k) {
        ETag[k]=0;
        ENumGmsh[k] = 0;
      }

      /// save info on triangular/tetrahedral elements
      int kTetra = 0;
      for (int k = 0; k < allElements; ++k) {
        meshfile >> eGmshNum >> eGeoType >> dummy;
        if (eGeoType==4){
          ENumGmsh[kTetra] = eGmshNum;           /// save the gmsh number
          meshfile >> ETag[kTetra] >> dummy;     /// save the element group (physical gmsh label)
          for(int v=0; v<mesh->Nverts; ++v)
            meshfile >> mesh->EToV[kTetra][v];         /// save the associated vertices
          kTetra++;
        }
        getline(meshfile, line);
      }

      int sk = 0, v;
      for (int k = 0; k < mesh->K; ++k){

        /* correct to 0-index */
        --(mesh->EToV[k][0]);
        --(mesh->EToV[k][1]);
        --(mesh->EToV[k][2]);
        --(mesh->EToV[k][3]);

        for(v=0;v<mesh->Nverts;++v){
          mesh->GX(k,v) = VX[mesh->EToV[k][v]];
          mesh->GY(k,v) = VY[mesh->EToV[k][v]];
	  mesh->GZ(k,v) = VZ[mesh->EToV[k][v]];
        }
	
      }

      std::cout << "BBDG: " << allVertices << " vertices and " << allElements << " elements in GMSH file " << std::endl;
      std::cout << "      (" << KTri << " triangles and " << KTetra << " tetrahedron)" << std::endl;
      std::cout << "mesh->K = " << mesh->K << std::endl;

    }
  }

#if 1
  /// correct orientation of negative elements by permuting their vertices
  for(int k=0; k<mesh->K; k++){

    double x1 = VX[mesh->EToV[k][0]], x2 = VX[mesh->EToV[k][1]], x3 = VX[mesh->EToV[k][2]], x4 = VX[mesh->EToV[k][3]];
    double y1 = VY[mesh->EToV[k][0]], y2 = VY[mesh->EToV[k][1]], y3 = VY[mesh->EToV[k][2]], y4 = VY[mesh->EToV[k][3]];
    double z1 = VZ[mesh->EToV[k][0]], z2 = VZ[mesh->EToV[k][1]], z3 = VZ[mesh->EToV[k][2]], z4 = VZ[mesh->EToV[k][3]];

    double xr = 0.5*(x2-x1), yr = 0.5*(y2-y1), zr = 0.5*(z2-z1);
    double xs = 0.5*(x3-x1), ys = 0.5*(y3-y1), zs = 0.5*(z3-z1);
    double xt = 0.5*(x4-x1), yt = 0.5*(y4-y1), zt = 0.5*(z4-z1);

    double J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

    if(J<0){
      //// swap two vertices to force positive Jacobian
      int tmp = mesh->EToV[k][3];
      mesh->EToV[k][3] = mesh->EToV[k][2];
      mesh->EToV[k][2] = tmp;
    }
  }
#endif

#define WRITE_MATLAB 0 // optional write to matlab
#if WRITE_MATLAB
  std::ofstream myfile;
  myfile.open("matlabMesh.m");
  myfile << "EToV = zeros(" << mesh->K << "," << p_Nfaces << ");" << std::endl;
  for(int e=0; e<mesh->K; ++e){
    myfile << "EToV(" << e+1 << ",:) = [";
    for(int f = 0; f< p_Nfaces; ++f){
      myfile << mesh->EToV[e][f]+1 << " ";
    }
    myfile << "];" << std::endl;
  }
  //  myfile << "VX = zeros(1," << mesh->Nv << ");" << std::endl;
  myfile << "VX = [";
  for(int i = 0; i < mesh->Nv; ++i){
    myfile << VX[i] << " ";
  }
  myfile << "];" << std::endl;
  myfile << "VY = [";
  for(int i = 0; i < mesh->Nv; ++i){
    myfile << VY[i] << " ";
  }
  myfile << "];" << std::endl;
  myfile << "VZ = [";
  for(int i = 0; i < mesh->Nv; ++i){
    myfile << VZ[i] << " ";
  }
  myfile << "];" << std::endl;

  myfile << "K = " << mesh->K << ";"<< std::endl;

  myfile.close();

#endif

  return mesh;

}

void PrintMesh ( Mesh *mesh ){
  int n;
  printf("Mesh data: \n");
  printf("\n K = %d\n", mesh->K);
  printf("\n Nv = %d\n", mesh->Nv);
  printf("\n Nverts = %d\n", mesh->Nverts);
  printf("\n Node coordinates = \n");
  printf("\n Element to vertex connectivity = \n");
  for(n=0;n<mesh->K;++n){
    printf("%d: %d %d %d %d\n", n,
	   mesh->EToV[n][0], mesh->EToV[n][1],
	   mesh->EToV[n][2], mesh->EToV[n][3]);
  }

}

void GeometricFactors3d(Mesh *mesh, int k,
                        double *drdx, double *dsdx, double *dtdx,
                        double *drdy, double *dsdy, double *dtdy,
                        double *drdz, double *dsdz, double *dtdz,
                        double *J){

  double x1 = mesh->GX(k,0), y1 =  mesh->GY(k,0), z1 =  mesh->GZ(k,0);
  double x2 = mesh->GX(k,1), y2 =  mesh->GY(k,1), z2 =  mesh->GZ(k,1);
  double x3 = mesh->GX(k,2), y3 =  mesh->GY(k,2), z3 =  mesh->GZ(k,2);
  double x4 = mesh->GX(k,3), y4 =  mesh->GY(k,3), z4 =  mesh->GZ(k,3);

  /* compute geometric factors of the following afine map */
  /* x = 0.5*( (-1-r-s-t)*x1 + (1+r)*x2 + (1+s)*x3 + (1+t)*x4) */
  /* y = 0.5*( (-1-r-s-t)*y1 + (1+r)*y2 + (1+s)*y3 + (1+t)*y4) */
  /* z = 0.5*( (-1-r-s-t)*z1 + (1+r)*z2 + (1+s)*z3 + (1+t)*z4) */

  double dxdr = (x2-x1)/2,  dxds = (x3-x1)/2, dxdt = (x4-x1)/2;
  double dydr = (y2-y1)/2,  dyds = (y3-y1)/2, dydt = (y4-y1)/2;
  double dzdr = (z2-z1)/2,  dzds = (z3-z1)/2, dzdt = (z4-z1)/2;

  *J =
     dxdr*(dyds*dzdt-dzds*dydt)
    -dydr*(dxds*dzdt-dzds*dxdt)
    +dzdr*(dxds*dydt-dyds*dxdt);

  *drdx  =  (dyds*dzdt - dzds*dydt)/(*J);
  *drdy  = -(dxds*dzdt - dzds*dxdt)/(*J);
  *drdz  =  (dxds*dydt - dyds*dxdt)/(*J);

  *dsdx  = -(dydr*dzdt - dzdr*dydt)/(*J);
  *dsdy  =  (dxdr*dzdt - dzdr*dxdt)/(*J);
  *dsdz  = -(dxdr*dydt - dydr*dxdt)/(*J);

  *dtdx  =  (dydr*dzds - dzdr*dyds)/(*J);
  *dtdy  = -(dxdr*dzds - dzdr*dxds)/(*J);
  *dtdz  =  (dxdr*dyds - dydr*dxds)/(*J);

  //  if(*J<1e-10)
    //    printf("warning: J = %lg on element %d\n", *J,k);

}

void Normals3d(Mesh *mesh, int k,
	       double *nx, double *ny, double *nz, double *sJ){

  int f;

  double drdx, dsdx, dtdx;
  double drdy, dsdy, dtdy;
  double drdz, dsdz, dtdz;
  double J;

  GeometricFactors3d(mesh, k,
		     &drdx, &dsdx, &dtdx,
		     &drdy, &dsdy, &dtdy,
		     &drdz, &dsdz, &dtdz,
		     &J);

  nx[0] = -dtdx; nx[1] = -dsdx; nx[2] =  drdx + dsdx + dtdx; nx[3] = -drdx;
  ny[0] = -dtdy; ny[1] = -dsdy; ny[2] =  drdy + dsdy + dtdy; ny[3] = -drdy;
  nz[0] = -dtdz; nz[1] = -dsdz; nz[2] =  drdz + dsdz + dtdz; nz[3] = -drdz;

  for(f=0;f<4;++f){
    sJ[f] = sqrt(nx[f]*nx[f]+ny[f]*ny[f]+nz[f]*nz[f]);
    nx[f] /= sJ[f];
    ny[f] /= sJ[f];
    nz[f] /= sJ[f];
    sJ[f] *= J;
  }
}

typedef struct {
  int k1, f1;
  int va, vb, vc;
}face3d;

int compare_pairs3d(const void *obj1, const void *obj2){

  face3d *e1 = (face3d*) obj1;
  face3d *e2 = (face3d*) obj2;

  int a1 = e1->va, b1 = e1->vb, c1 = e1->vc;
  int a2 = e2->va, b2 = e2->vb, c2 = e2->vc;

  int va1, vb1, vc1, va2, vb2, vc2;

  va1 = min(a1, min(b1, c1));
  vc1 = max(a1, max(b1, c1));

  if(va1!=a1 && vc1!=a1) vb1=a1;
  if(va1!=b1 && vc1!=b1) vb1=b1;
  if(va1!=c1 && vc1!=c1) vb1=c1;

  va2 = min(a2, min(b2, c2));
  vc2 = max(a2, max(b2, c2));

  if(va2!=a2 && vc2!=a2) vb2=a2;
  if(va2!=b2 && vc2!=b2) vb2=b2;
  if(va2!=c2 && vc2!=c2) vb2=c2;

  if(vc1<vc2)
    return -1;
  else if(vc1>vc2)
    return 1;
  else if(vb1<vb2)
    return -1;
  else if(vb1>vb2)
    return 1;
  else if(va1<va2)
    return -1;
  else if(va1>va2)
    return 1;

  return 0;

}

void FacePair3d(Mesh *mesh){

  int K = mesh->K;
  int Nfaces = mesh->Nfaces;
  int Nverts = mesh->Nverts;

  int **EToV = mesh->EToV;

  const int vnum[4][3] = { {0,1,2}, {0,1,3}, {1,2,3}, {0,2,3} };

  int n, k, e, sk, v;

  // tag faces by their sorted vertices
  face3d *myfaces = (face3d*) calloc(K*Nfaces, sizeof(face3d));

  sk = 0;
  for(k=0;k<K;++k){
    for(e=0;e<Nfaces;++e){
      int a1 = EToV[k][vnum[e][0]];
      int b1 = EToV[k][vnum[e][1]];
      int c1 = EToV[k][vnum[e][2]];

      myfaces[sk].k1 = k; myfaces[sk].f1 = e;

      int va1, vb1, vc1, va2, vb2, vc2;

      va1 = min(a1, min(b1, c1));
      vc1 = max(a1, max(b1, c1));

      if(va1!=a1 && vc1!=a1) vb1=a1;
      if(va1!=b1 && vc1!=b1) vb1=b1;
      if(va1!=c1 && vc1!=c1) vb1=c1;

      myfaces[sk].va = va1;
      myfaces[sk].vb = vb1;
      myfaces[sk].vc = vc1;
      ++sk;
    }
  }

  qsort(myfaces, K*Nfaces, sizeof(face3d), compare_pairs3d);

  mesh->EToE = BuildIntMatrix(K, Nfaces);
  mesh->EToF = BuildIntMatrix(K, Nfaces);

  int id, k1, k2, f1, f2, p1, p2;
  sk = 0;

  for(n=0; n < K*Nfaces-1; ++n){

    int e1  = myfaces[n].k1;   int e2  = myfaces[n+1].k1;
    int f1  = myfaces[n].f1;   int f2  = myfaces[n+1].f1;
    int va1 = myfaces[n].va;   int va2 = myfaces[n+1].va;
    int vb1 = myfaces[n].vb;   int vb2 = myfaces[n+1].vb;
    int vc1 = myfaces[n].vc;   int vc2 = myfaces[n+1].vc;

    // if faces match
    if (va1==va2 && vb1==vb2 && vb1==vb2){

      mesh->EToE[e1][f1] = e2;
      mesh->EToF[e1][f1] = f2;
      mesh->EToE[e2][f2] = e1;
      mesh->EToF[e2][f2] = f1;
      ++n; // skip to next face

    }else{ // otherwise, set as boundary face

      mesh->EToE[e1][f1] = e1;
      mesh->EToF[e1][f1] = f1;

    }
  }


  free(myfaces);
}
