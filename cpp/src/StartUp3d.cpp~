#include "fem.h"
#include "Basis.h"

void StartUp3d(Mesh *mesh){

  // default to planar elements at startup
  mesh->KCurved = 0;
  mesh->KPlanar = mesh->K;

  VectorXi KlistCurved(1); KlistCurved.fill(0);
  mesh->KlistCurved = KlistCurved;
  VectorXi KlistCurvedPlusNbrs(1); KlistCurvedPlusNbrs.fill(0);
  mesh->KlistCurvedPlusNbrs = KlistCurvedPlusNbrs;

  VectorXi KlistPlanar(mesh->K);
  for (int e = 0; e < mesh->K; ++e){
    KlistPlanar(e) = e;
  }
  mesh->KlistPlanar = KlistPlanar;

  int n, m, k;

  // ======= computing points, dmats =======

  VectorXd r,s,t;
  Nodes3D(p_N,r,s,t);

  MatrixXd V = Vandermonde3D(p_N,r,s,t);
  mesh->V = V;
  MatrixXd Vr,Vs,Vt;
  GradVandermonde3D(p_N,r,s,t,Vr,Vs,Vt);
  MatrixXd Dr = mrdivide(Vr,V);
  MatrixXd Ds = mrdivide(Vs,V);
  MatrixXd Dt = mrdivide(Vt,V);

  //  cout << "r = " << endl << r << endl;
  //  cout << "s = " << endl << s << endl;
  //  cout << "t = " << endl << t << endl;

  // set up Fmask for face nodes
  MatrixXi Fmask(p_Nfp,p_Nfaces);
  int sk0 = 0, sk1 = 0, sk2 = 0, sk3 = 0;
  for(int i = 0; i < p_Np; ++i){
    if (fabs(1. + t(i)) < NODETOL){
      Fmask(sk0,0) = i; ++sk0;
    }
    if (fabs(1. + s(i)) < NODETOL){
      Fmask(sk1,1) = i; ++sk1;
    }
    if (fabs(1. + r(i) + s(i) + t(i)) < NODETOL){
      Fmask(sk2,2) = i; ++sk2;
    }
    if (fabs(1. + r(i)) < NODETOL){
      Fmask(sk3,3) = i; ++sk3;
    }
  }
  mesh->Fmask = Fmask;

  //cout << "Fmask = " << endl << Fmask << endl;
  //  cout << "p_Nfp = " << p_Nfp <<endl;

  // get face nodes
  MatrixXd rf,sf;
  rf.resize(p_Nfp,p_Nfaces);
  sf.resize(p_Nfp,p_Nfaces);
  for (int f = 0; f < p_Nfaces; ++f){
    for (int i = 0; i < p_Nfp; ++i){
      if (f==0){
	rf(i,f) = r(Fmask(i,f));  sf(i,f) = s(Fmask(i,f));
      }else if (f==1){
	rf(i,f) = r(Fmask(i,f));  sf(i,f) = t(Fmask(i,f));
      }else if (f==2){
	rf(i,f) = s(Fmask(i,f));  sf(i,f) = t(Fmask(i,f));
      }else if (f==3){
	rf(i,f) = s(Fmask(i,f));  sf(i,f) = t(Fmask(i,f));
      }
    }
  }

  int qN = min(21,2*p_N+1); // max order = 21
  int qNf = min(21,2*p_N+1);
  //int qN = min(21,4*p_N-3); // max order = 21
  //int qNf = min(21,4*p_N-2);

  VectorXd rq,sq,tq,wq;
  tet_cubature(qN,rq,sq,tq,wq); // can increase quadrature degree up to 21
  mesh->rq = rq;
  mesh->sq = sq;
  mesh->tq = tq;
  mesh->wq = wq;
  mesh->Nq = wq.rows();

  VectorXd rqtri,sqtri,wqtri;
  tri_cubature(qNf,rqtri,sqtri,wqtri);
  mesh->Nfq = wqtri.rows();
  mesh->wfqFace = wqtri;

  printf("qN = %d, qNf = %d\n",qN,qNf);

  MatrixXd Vfqtmp = Vandermonde2D(p_N,rqtri,sqtri);
  MatrixXd Vtri = Vandermonde2D(p_N,rf.col(0),sf.col(0)); // use bottom face of tet
  MatrixXd VfqFace = mrdivide(Vfqtmp,Vtri);
  mesh->VfqFace = VfqFace; // interpolate from tri nodes to triq nodes
  //  cout << "VfqFace = " << endl << VfqFace << endl;
  //  printf("VfqFace dimensions = %d by %d\n",VfqFace.rows(),VfqFace.cols());
  //  printf("Nfq = %d, NfqNfaces = %d\n",mesh->Nfq, mesh->Nfq*p_Nfaces);

  // all faces
  VectorXd rfq,sfq,tfq,wfq;
  int Nfq = mesh->Nfq;
  rfq.resize(Nfq*p_Nfaces);
  sfq.resize(Nfq*p_Nfaces);
  tfq.resize(Nfq*p_Nfaces);
  wfq.resize(Nfq*p_Nfaces);

  int foff = 0;
  // face 1
  rfq.middleRows(foff,Nfq) = rqtri;
  sfq.middleRows(foff,Nfq) = sqtri;
  tfq.middleRows(foff,Nfq).fill(-1.0);
  wfq.middleRows(foff,Nfq) = wqtri;
  foff += Nfq;

  // face 2
  rfq.middleRows(foff,Nfq) = rqtri;
  sfq.middleRows(foff,Nfq).fill(-1.0);
  tfq.middleRows(foff,Nfq) = sqtri;
  wfq.middleRows(foff,Nfq) = wqtri;
  foff += Nfq;

  // face 3
  rfq.middleRows(foff,Nfq) = -(1+rqtri.array()+sqtri.array());
  sfq.middleRows(foff,Nfq) = rqtri;
  tfq.middleRows(foff,Nfq) = sqtri;
  wfq.middleRows(foff,Nfq) = wqtri;
  foff += Nfq;

  // face 4
  rfq.middleRows(foff,Nfq).fill(-1);
  sfq.middleRows(foff,Nfq) = rqtri;
  tfq.middleRows(foff,Nfq) = sqtri;
  wfq.middleRows(foff,Nfq) = wqtri;

  /*
  cout << "Fmask = [" << endl << Fmask << "];" << endl;
  cout << "rfq = [" << endl << rfq << "];" << endl;
  cout << "sfq = [" << endl << sfq << "];" << endl;
  cout << "tfq = [" << endl << tfq << "];" << endl;
  */

  Vfqtmp = Vandermonde3D(p_N,rfq,sfq,tfq);
  MatrixXd Vfq = mrdivide(Vfqtmp,V);
  mesh->rfq = rfq;
  mesh->sfq = sfq;
  mesh->tfq = tfq;
  mesh->Vfq = Vfq; // interp from nodes to surface cubature nodes
  mesh->wfq = wfq;
  // save nodes
#if 0
  mesh->r = BuildVector(p_Np);
  mesh->s = BuildVector(p_Np);
  mesh->t = BuildVector(p_Np);
  for (int i = 0; i < p_Np; ++i){
    mesh->r[i] = r(i);
    mesh->s[i] = s(i);
    mesh->t[i] = t(i);
  }
#else
  mesh->r = r;
  mesh->s = s;
  mesh->t = t;
#endif

  // save Fmask as a C array
  mesh->FmaskC = BuildIntMatrix(p_Nfaces, p_Nfp);
  for(int n=0;n<p_Nfaces;++n){
    for(int m=0;m<p_Nfp;++m){
      mesh->FmaskC[n][m] = Fmask(m,n);
    }
  }

  // save quadrature nodes (need weights to compute errors)
  mesh->Dr = Dr;
  mesh->Ds = Ds;
  mesh->Dt = Dt;

  // build lift matrix
  MatrixXd Emat(p_Np,p_Nfp*p_Nfaces);
  Emat.fill(0.0);
  for (int f = 0; f < p_Nfaces; ++f){
    MatrixXd Vface = Vandermonde2D(p_N,rf.col(f),sf.col(f));
    MatrixXd massFace = (Vface*Vface.transpose()).inverse();

    for (int i = 0; i < p_Nfp; ++i){
      for (int j = 0; j < p_Nfp; ++j){
	int idr = Fmask(i,f);
	int idc = j + f*p_Nfp;
	Emat(idr,idc) += massFace(i,j);
      }
    }
  }
  MatrixXd LIFT = (V*V.transpose())*Emat; // V*V' = inv(M) for nodal
#if 0
  mesh->LIFT = BuildMatrix(p_Np, p_Nfp*p_Nfaces);
  for(n=0;n<p_Np;++n){
    for(m=0;m<p_Nfp*p_Nfaces;++m){
      mesh->LIFT[n][m] = LIFT(n,m);
    }
  }
#else
  mesh->LIFT = LIFT;
#endif

  // quadrature matrix
  MatrixXd Vqtmp = Vandermonde3D(p_N,rq,sq,tq);
  MatrixXd Vq = mrdivide(Vqtmp,V);
#if 0
  mesh->Vq = BuildMatrix(mesh->Nq,p_Np);
  for(int i = 0; i < mesh->Nq; ++i){
    for (int j = 0; j < p_Np; ++j){
      mesh->Vq[i][j] = Vq(i,j);
    }
  }
#else
  mesh->Vq = Vq;
#endif

  // Bernstein conversion to nodal matrix
  MatrixXd VB = BernTet(p_N,r,s,t);
  MatrixXd invVB = VB.inverse();
#if 0
  // save Bernstein VDMs [JC]
  mesh->VB = BuildMatrix(p_Np, p_Np);
  mesh->invVB = BuildMatrix(p_Np, p_Np);
  for(n=0;n<p_Np;++n){
    for(m=0;m<p_Np;++m){
      mesh->VB[n][m] = VB(n,m);
      mesh->invVB[n][m] = invVB(n,m);
    }
  }
#else
  mesh->VB = VB;
  mesh->invVB = invVB;
#endif

  // Bernstein deriv matrices
  MatrixXd VB1,VB2,VB3,VB4;
  GradBernTet(p_N,r,s,t,VB1,VB2,VB3,VB4);
  MatrixXd DB1 = mldivide(VB,VB1);
  MatrixXd DB2 = mldivide(VB,VB2);
  MatrixXd DB3 = mldivide(VB,VB3);
  MatrixXd DB4 = mldivide(VB,VB4);

  MatrixXi D1_ids,D2_ids,D3_ids,D4_ids;
  MatrixXd D_vals;
  get_sparse_ids(DB1,D1_ids,D_vals);
  get_sparse_ids(DB2,D2_ids,D_vals);
  get_sparse_ids(DB3,D3_ids,D_vals);
  get_sparse_ids(DB4,D4_ids,D_vals);

  // bernstein sparse D ops
  mesh->D1_ids = BuildIntMatrix(p_Np,4);
  mesh->D2_ids = BuildIntMatrix(p_Np,4);
  mesh->D3_ids = BuildIntMatrix(p_Np,4);
  mesh->D4_ids = BuildIntMatrix(p_Np,4);
  mesh->D_vals = BuildMatrix(p_Np,4);
  for(int i = 0; i < p_Np; ++i){
    for(int j = 0; j < 4; ++j){
      mesh->D1_ids[i][j] = 0;
      mesh->D2_ids[i][j] = 0;
      mesh->D3_ids[i][j] = 0;
      mesh->D4_ids[i][j] = 0;
      mesh->D_vals[i][j] = 0.0;
      if (j < D1_ids.cols()){
	mesh->D1_ids[i][j] = D1_ids(i,j);
	mesh->D2_ids[i][j] = D2_ids(i,j);
	mesh->D3_ids[i][j] = D3_ids(i,j);
	mesh->D4_ids[i][j] = D4_ids(i,j);
	mesh->D_vals[i][j] = D_vals(i,j);
      }
    }
  }

  // bernstein lift operators
  MatrixXd Vface = BernTri(p_N,rqtri,sqtri);
  MatrixXd massFace = Vface.transpose() * mesh->wfqFace.asDiagonal() * Vface; // mass matrix by quadrature

  //cout << "rfq = " << endl << rfq << endl;
  //cout << "sfq = " << endl << sfq << endl;
  //cout << "bern vface = " << endl << Vface << endl;
  //cout << "bern massface = " << endl << massFace << endl;

  MatrixXd EmatB(p_Np,p_Nfp*p_Nfaces);
  EmatB.fill(0.0);
  for (int f = 0; f < p_Nfaces; ++f){

    for (int i = 0; i < p_Nfp; ++i){
      for (int j = 0; j < p_Nfp; ++j){
	int idr = Fmask(i,f);
	int idc = j + f*p_Nfp;
	EmatB(idr,idc) += massFace(i,j);
      }
    }
  }
  MatrixXd VBq = BernTet(p_N,rq,sq,tq);
  MatrixXd MB = VBq.transpose() * wq.asDiagonal() * VBq;
  MatrixXd LIFTB = mldivide(MB,EmatB);
  // cout << "MB = " << endl << MB << endl;
  // cout << "EmatB = " << endl << EmatB << endl;
  // cout << "LIFTB = " << endl << LIFTB << endl;

  mesh->cEL.resize(p_N+1);
  MatrixXd EEL1(p_Np, p_Nfp);
  MatrixXd VBTri = BernTri(p_N,rqtri,sqtri);
  int off = 0;
  for (int i = 0; i <= p_N; ++i){
    int NpTri = (p_N+1)*(p_N+2)/2;
    int NpElev = (p_N-i+1)*(p_N-i+2)/2;
    MatrixXd Ei(NpTri,NpElev);
    MatrixXd VBi = BernTri(p_N-i,rqtri,sqtri);
    Ei = mldivide(VBTri,VBi);
    double cEL = pow(-1,i)*nchoosek(p_N,i)/(1+i);
    mesh->cEL(i) = cEL;
    EEL1.middleRows(off,NpElev) = cEL*Ei.transpose();
    off += NpElev;
  }
  //cout << "EEL1 = " << endl << EEL1 << endl;

  // get permutation of other faces
  MatrixXi rowperm(p_Np,p_Nfaces-1);
  off = p_Nfp;
  MatrixXd LIFTBf1 = LIFTB.leftCols(p_Nfp);
  double maxVal = LIFTB.array().abs().maxCoeff();

  for (int f = 1; f < p_Nfaces; ++f){
    MatrixXd LIFTBf = LIFTB.middleCols(off,p_Nfp);
    for (int i = 0; i < p_Np; ++i){
      int match_found = 0;
      int row = 0;
      while(match_found==0){
	ArrayXd diff = (LIFTBf1.row(row).array() - LIFTBf.row(i).array()).abs();

	if (diff.maxCoeff() < NODETOL*maxVal){ // WARNING: BB can be sensitive to this tolerance!!
	  match_found = 1;
	}else{
	  ++row;
	}
      }
      rowperm(i,f-1) = row;
    }
    off += p_Nfp;
  }
  //cout << "rowperm = " << endl << rowperm << endl;
  mesh->faceVolPerm = rowperm;

  // make full lift reduction
  MatrixXd EEL(p_Np,p_Nfaces*p_Nfp);
  EEL.leftCols(p_Nfp) = EEL1;
  off = p_Nfp;
  for (int f = 1; f < p_Nfaces;++f){
    MatrixXd EELperm(p_Np,p_Nfp);
    for (int i = 0; i < p_Np; ++i){
      EELperm.row(i) = EEL1.row(rowperm(i,f-1));
    }
    EEL.middleCols(off,p_Nfp) = EELperm;
    off += p_Nfp;
  }

  //cout << "EEL = " << endl << EEL << endl;
  MatrixXi EEL_ids;
  MatrixXd EEL_vals;
  get_sparse_ids(EEL,EEL_ids,EEL_vals);

  // decomposition of LIFT operator
  mesh->EEL_nnz = p_Nfp + 3; // max nonzeros per row
  mesh->EEL_ids = BuildIntMatrix(p_Np,mesh->EEL_nnz);
  mesh->EEL_vals = BuildMatrix(p_Np,mesh->EEL_nnz);
  for(int i = 0; i < p_Np; ++i){
    for (int j = 0; j < mesh->EEL_nnz; ++j){
      mesh->EEL_vals[i][j] = 0.0;
      mesh->EEL_ids[i][j] = 0;
      if (j < EEL_ids.cols()){
	mesh->EEL_vals[i][j] = EEL_vals(i,j);
	mesh->EEL_ids[i][j] = EEL_ids(i,j);
      }
    }
  }

  // Nfp x Nfp L0 operator
  MatrixXd L0 = LIFTB.topLeftCorner(p_Nfp,p_Nfp);
  MatrixXi L0_ids;
  MatrixXd L0_vals;
  get_sparse_ids(L0,L0_ids,L0_vals);

  int L0_nnz = min(p_Nfp,7); // fixed L0 nnz per row
  mesh->L0_ids = BuildIntMatrix(p_Nfp,L0_nnz);
  mesh->L0_vals = BuildMatrix(p_Nfp,L0_nnz);
  for(int i = 0; i < p_Nfp; ++i){
    for (int j = 0; j < L0_nnz; ++j){
      mesh->L0_vals[i][j] = 0.0;
      mesh->L0_ids[i][j]  = 0;
      if (j < L0_ids.cols()){
	mesh->L0_vals[i][j] = L0_vals(i,j);
	mesh->L0_ids[i][j]  = L0_ids(i,j);
      }
    }
  }

  // store face reduction sparsely - col indices + vals
  MatrixXi vol_ids(p_Np,p_Nfaces);
  for (int i = 0; i < p_Np; ++i){
    vol_ids(i,0) = i;
    vol_ids(i,1) = rowperm(i,0);
    vol_ids(i,2) = rowperm(i,1);
    vol_ids(i,3) = rowperm(i,2);
  }
  mesh->vol_ids = vol_ids;

  MatrixXi slice_ids(p_Np,4); // assumes nfaces = 4

#define MAPJ(NN,JJ)                             \
  ((NN+1)*JJ - JJ*(JJ+1)/2 + JJ)
#define MAPK(NN,KK)                                             \
  (KK*(3*NN*(NN - KK + 4) + KK*(KK - 6) + 11)/6)

  int sk = 0;
  for (int slice = 0; slice <= p_N; ++slice){
    for (int j = 0; j <= p_N-slice; ++j){
      for (int i = 0; i <= p_N-slice-j; ++i){
	int l = p_N-i-j-k;

	// get vol_id
	// face 1: k = slice
	int slice_offset = MAPK(p_N,slice);
	int line_offset = MAPJ(p_N-slice,j);
	slice_ids(sk,0) = i + line_offset + slice_offset;

	// face 2: j = slice, set k = j
	slice_offset = MAPK(p_N,j);
	line_offset = MAPJ(p_N-j,slice);
	slice_ids(sk,1) = i + line_offset + slice_offset;

	// face 3: set l = slice, k = j, j = i, i = N-j-k-l
	slice_offset = MAPK(p_N,j);
	line_offset  = MAPJ(p_N-j,i);
	slice_ids(sk,2) = (p_N-slice-i-j) + line_offset + slice_offset;

	// face 4: i = slice, k = j, j = i
	slice_offset = MAPK(p_N,j);
	line_offset = MAPJ(p_N-j,i);
	slice_ids(sk,3) = slice + line_offset + slice_offset;

	++sk;
      }
    }
  }
  mesh->slice_ids = slice_ids;
  //cout << "slice_ids = " << endl << slice_ids << endl;

  vector<MatrixXd> EiVec;
  for (int i = 0; i < p_N; ++i){
    MatrixXd VB1 = BernTri(p_N-i,rqtri,sqtri);
    MatrixXd VB2 = BernTri(p_N-i-1,rqtri,sqtri);
    MatrixXd Ei = mldivide(VB1,VB2);
    EiVec.push_back(Ei); // store for later
    //cout << "Ei = " << endl << Ei << endl;
  }

  VectorXd EEL_val_vec;
  VectorXi EEL_id_vec;
  //cout << "EEL has " << EEL_val_vec.rows() << " rows" << endl;
  for (int i = 0; i < p_N; ++i){
    MatrixXd EiTr = EiVec[i].transpose();
    MatrixXi Ei_ids;
    MatrixXd Ei_vals;
    get_sparse_ids(EiTr,Ei_ids,Ei_vals);

    // expand to () x 3 arrays
    if (Ei_vals.cols() < 3){
      MatrixXd zeroCols(Ei_vals.rows(),3-Ei_vals.cols());
      zeroCols.fill(0.0);
      MatrixXi zeroIds(Ei_ids.rows(),3-Ei_ids.cols());
      zeroIds.fill(0);
      MatrixXd joinedVals(Ei_vals.rows(),3);
      joinedVals << Ei_vals, zeroCols;
      MatrixXi joinedIds(Ei_ids.rows(),3);
      joinedIds << Ei_ids, zeroIds;

      Ei_vals = joinedVals;
      Ei_ids = joinedIds;
    }
    //cout << "Ei_vals = " << endl << Ei_vals << endl;
    //cout << "Ei_ids = " << endl << Ei_ids << endl;

    // flatten matrix to vector
    VectorXd Ei_vals_vec(Map<VectorXd>(Ei_vals.data(), Ei_vals.cols()*Ei_vals.rows()));
    VectorXi Ei_ids_vec(Map<VectorXi>(Ei_ids.data(), Ei_ids.cols()*Ei_ids.rows()));

    VectorXd joined(EEL_val_vec.rows() + Ei_vals_vec.rows());
    joined << EEL_val_vec, Ei_vals_vec;
    EEL_val_vec = joined;

    VectorXi id_joined(EEL_id_vec.rows() + Ei_ids_vec.rows());
    id_joined << EEL_id_vec, Ei_ids_vec;
    EEL_id_vec = id_joined;

  }
  mesh->EEL_val_vec = EEL_val_vec;
  mesh->EEL_id_vec = EEL_id_vec;
  //cout << "EEL val vec = " << endl << EEL_val_vec << endl;
  //cout << "EEL id vec = " << endl << EEL_id_vec << endl;


  // ====================================================

  // low storage RK coefficients
  mesh->rk4a = BuildVector(5);
  mesh->rk4a[0] =              0.0;
  mesh->rk4a[1] =  -567301805773.0 / 1357537059087.0;
  mesh->rk4a[2] = -2404267990393.0 / 2016746695238.0;
  mesh->rk4a[3] = -3550918686646.0 / 2091501179385.0;
  mesh->rk4a[4] = -1275806237668.0 /  842570457699.0;

  mesh->rk4b = BuildVector(5);
  mesh->rk4b[0] =  1432997174477.0 /  9575080441755.0;
  mesh->rk4b[1] =  5161836677717.0 / 13612068292357.0;
  mesh->rk4b[2] =  1720146321549.0 /  2090206949498.0;
  mesh->rk4b[3] =  3134564353537.0 /  4481467310338.0;
  mesh->rk4b[4] =  2277821191437.0 / 14882151754819.0;

  mesh->rk4c = BuildVector(6);
  mesh->rk4c[0] =              0.0;
  mesh->rk4c[1] =  1432997174477.0 / 9575080441755.0;
  mesh->rk4c[2] =  2526269341429.0 / 6820363962896.0;
  mesh->rk4c[3] =  2006345519317.0 / 3224310063776.0;
  mesh->rk4c[4] =  2802321613138.0 / 2924317926251.0;
  mesh->rk4c[5] =              1.0;

  // map coordinates
  //mesh->x = BuildMatrix(mesh->K, p_Np);
  //mesh->y = BuildMatrix(mesh->K, p_Np);
  //mesh->z = BuildMatrix(mesh->K, p_Np);
  mesh->x.resize(p_Np,mesh->K);
  mesh->y.resize(p_Np,mesh->K);
  mesh->z.resize(p_Np,mesh->K);

  for(k=0;k<mesh->K;++k){
    for(n=0;n<p_Np;++n){
      //dfloat r = mesh->r[n];
      //dfloat s = mesh->s[n];
      //dfloat t = mesh->t[n];
      double r = mesh->r(n);
      double s = mesh->s(n);
      double t = mesh->t(n);

      mesh->x(n,k) = 0.5*( -mesh->GX(k,0)*(r+s+t+1) +
			   mesh->GX(k,1)*(1.+r) +
			   mesh->GX(k,2)*(1.+s) +
			   mesh->GX(k,3)*(1.+t)
			   );

      mesh->y(n,k) = 0.5*( -mesh->GY(k,0)*(r+s+t+1) +
			   mesh->GY(k,1)*(1.+r) +
			   mesh->GY(k,2)*(1.+s) +
			   mesh->GY(k,3)*(1.+t)
			   );

      mesh->z(n,k) = 0.5*( -mesh->GZ(k,0)*(r+s+t+1) +
			   mesh->GZ(k,1)*(1.+r) +
			   mesh->GZ(k,2)*(1.+s) +
			   mesh->GZ(k,3)*(1.+t)
			    );
    }
  }

  // compute geofacs for straight-sided mesh [JC]
  mesh->rx.resize(p_Np,mesh->K);
  mesh->sx.resize(p_Np,mesh->K);
  mesh->tx.resize(p_Np,mesh->K);
  mesh->ry.resize(p_Np,mesh->K);
  mesh->sy.resize(p_Np,mesh->K);
  mesh->ty.resize(p_Np,mesh->K);
  mesh->rz.resize(p_Np,mesh->K);
  mesh->sz.resize(p_Np,mesh->K);
  mesh->tz.resize(p_Np,mesh->K);
  mesh->J.resize( p_Np,mesh->K);

  mesh->nx.resize(p_Nfp*p_Nfaces,mesh->K);
  mesh->ny.resize(p_Nfp*p_Nfaces,mesh->K);
  mesh->nz.resize(p_Nfp*p_Nfaces,mesh->K);
  mesh->sJ.resize(p_Nfp*p_Nfaces,mesh->K);

  double drdx, dsdx, dtdx;
  double drdy, dsdy, dtdy;
  double drdz, dsdz, dtdz, J;
  double *nxk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  double *nyk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  double *nzk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  double *sJk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(mesh->Nfaces);
  for(int k=0;k<mesh->K;++k){

    GeometricFactors3d(mesh, k,
		       &drdx, &dsdx, &dtdx,
		       &drdy, &dsdy, &dtdy,
		       &drdz, &dsdz, &dtdz, &J);

    mesh->rx.col(k).fill(drdx);
    mesh->sx.col(k).fill(dsdx);
    mesh->tx.col(k).fill(dtdx);
    mesh->ry.col(k).fill(drdy);
    mesh->sy.col(k).fill(dsdy);
    mesh->ty.col(k).fill(dtdy);
    mesh->rz.col(k).fill(drdz);
    mesh->sz.col(k).fill(dsdz);
    mesh->tz.col(k).fill(dtdz);
    mesh->J.col(k).fill(J);

    Normals3d(mesh, k, nxk, nyk, nzk, sJk);

    int off = 0;
    for(int f=0;f<mesh->Nfaces;++f){

      double Fscale = sJk[f]/J; //sJk[f]/(2.*J);
      double nx = nxk[f];
      double ny = nyk[f];
      double nz = nzk[f];
      for (int i = 0; i < p_Nfp; ++i){
	mesh->nx(i+off,k) = nx;
	mesh->ny(i+off,k) = ny;
	mesh->nz(i+off,k) = nz;
	mesh->sJ(i+off,k) = sJk[f];
      }
      off += p_Nfp;
    }
  }

  // build node-node connectivity maps
  void BuildMaps3d(Mesh *mesh);
  BuildMaps3d(mesh);

}

void BuildMaps3d(Mesh *mesh){

  //printf("Hello %d\n", 1002);

  int K       = mesh->K;
  int Nfaces  = mesh->Nfaces;

  mesh->vmapM = BuildIntVector(p_Nfp*p_Nfaces*K);
  mesh->vmapP = BuildIntVector(p_Nfp*p_Nfaces*K);

  mesh->mapM = BuildIntVector(p_Nfp*p_Nfaces*K);
  mesh->mapP = BuildIntVector(p_Nfp*p_Nfaces*K);

  int m;
  int k1,f1,p1,n1,id1, k2,f2,p2,n2,id2;

  double x1, y1,z1, x2, y2, z2, d12;

  double *nxk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(Nfaces);
  double *nyk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(Nfaces);
  double *nzk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(Nfaces);
  double *sJk = (double*) calloc(mesh->Nfaces,sizeof(double));//BuildVector(Nfaces);

  printf("Hello %d\n", 1001);

  /* first build local */
  for(k1=0;k1<K;++k1){

    /* get some information about the face geometries */
    Normals3d(mesh, k1, nxk, nyk, nzk, sJk);

    for(f1=0;f1<Nfaces;++f1){

      /* volume -> face nodes */
      for(n1=0;n1<p_Nfp;++n1){
	id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces;
	mesh->vmapM[id1] = mesh->FmaskC[f1][n1] + k1*p_Np;
	mesh->mapM[id1] = id1;
      }

      /* find neighbor */
      k2 = mesh->EToE[k1][f1];
      f2 = mesh->EToF[k1][f1];

      if(k1==k2){
	for(n1=0;n1<p_Nfp;++n1){
	  id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces;
	  mesh->vmapP[id1] = k1*p_Np + mesh->FmaskC[f1][n1];
	  mesh->mapP[id1] = mesh->mapM[id1];
	}
      }else{
	/* treat as boundary for the moment  */

	for(n1=0;n1<p_Nfp;++n1){
	  id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces;

	  x1 = mesh->x(mesh->FmaskC[f1][n1],k1);
	  y1 = mesh->y(mesh->FmaskC[f1][n1],k1);
	  z1 = mesh->z(mesh->FmaskC[f1][n1],k1);

	  for(n2=0;n2<p_Nfp;++n2){

	    id2 = n2+f2*p_Nfp+k2*p_Nfp*p_Nfaces;

	    x2 = mesh->x(mesh->FmaskC[f2][n2],k2);
	    y2 = mesh->y(mesh->FmaskC[f2][n2],k2);
	    z2 = mesh->z(mesh->FmaskC[f2][n2],k2);

	    /* find normalized distance between these nodes */
	    /* [ use sJk as a measure of edge length (ignore factor of 2) ] */
	    d12 = ((x1-x2)*(x1-x2) +
		   (y1-y2)*(y1-y2) +
		   (z1-z2)*(z1-z2)); /* /(sJk[f1]*sJk[f1]);  */
	    if(d12<NODETOL){
	      mesh->vmapP[id1] = k2*p_Np + mesh->FmaskC[f2][n2];
	      mesh->mapP[id1] = id2;
	      break;
	    }
	  }
	  if(n2==p_Nfp){
	    printf("LOST NODE on elem %d, face %d!!!\n",k1,f1);
	  }
	}
      }
    }
  }

#if 0
  for (int i = 0; i < p_Nfp*p_Nfaces;++i){
    printf("mesh->mapP[%d] = ",i);
    for (int e= 0; e < mesh->K; ++e){
      int id = i + e*p_Nfp*p_Nfaces;
      printf("%d ", mesh->mapP[id]);
    }
    printf("\n");
  }
#endif
  return;
}


// general: input nodes (ex: quadrature nodes), get map back
void BuildFaceNodeMaps(Mesh *mesh, MatrixXd xf, MatrixXd yf, MatrixXd zf,
		       MatrixXi &mapP){

  int K       = mesh->K;
  int Nfaces  = mesh->Nfaces;
  //  printf("In buildfacenodemaps Nfaces = %d\n",Nfaces);

  mapP.resize(xf.rows(),K);
  mapP.fill(-1);

  int Nfpts = xf.rows() / Nfaces; // assume same # qpts per face

  int m;
  int k1,p1,n1,f1,k2,p2,n2,f2;

  double x1, y1, z1, x2, y2, z2, d12;

  printf("Hello %d\n", 1337);

  // first build local
  for(k1=0;k1<K;++k1){

    for(f1=0;f1<Nfaces;++f1){

      // find neighbor
      k2 = mesh->EToE[k1][f1];
      f2 = mesh->EToF[k1][f1];

      if(k1==k2){
	for(int i = 0; i < Nfpts; ++i){
	  mapP(i + f1*Nfpts,k1) = i + f1*Nfpts + xf.rows()*k1;
	}
      }else{

        MatrixXd xM(Nfpts,3);
        MatrixXd xP(Nfpts,3);
        for(int i = 0; i < Nfpts; ++i){
	  int id1 = i + Nfpts * f1;
	  x1 = xf(id1,k1); y1 = yf(id1,k1); z1 = zf(id1,k1);
          int id2 = i + Nfpts * f2;
          x2 = xf(id2,k2); y2 = yf(id2,k2); z2 = zf(id2,k2);

          xM(i,0) = x1; xM(i,1) = y1; xM(i,2) = z1;
          xP(i,0) = x2; xP(i,1) = y2; xP(i,2) = z2;
        }

	for(int i = 0; i < Nfpts; ++i){

	  double minDist = 1000.0;

	  int id1 = i + Nfpts * f1;
	  x1 = xf(id1,k1); y1 = yf(id1,k1); z1 = zf(id1,k1);
          bool node_found = false;
	  for(int j = 0; j < Nfpts; ++j){

	    int id2 = j + Nfpts * f2;
	    x2 = xf(id2,k2); y2 = yf(id2,k2); z2 = zf(id2,k2);

	    // find distance between these nodes
	    d12 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
	    minDist = min(minDist,d12);
	    if (d12<NODETOL){
	      mapP(id1,k1) = id2 + xf.rows()*k2;
	      node_found = true;
	      break;
	    }
	  }
          if(!node_found){
            printf("BuildFaceNodeMaps: lost node %d on elems %d, %d. min dist = %g\n",i,k1,k2,minDist);
            cout << "xM = " << endl << xM << endl;
            cout << "xP = " << endl << xP << endl;
          }
	}
      }

    }// faces
  }// k


  return;

}


void InitQuadratureArrays(Mesh *mesh){

  int Nq = mesh->Nq;
  int Nfq = mesh->Nfq;
  int Nfaces = mesh->Nfaces;

  // matlab-esque: requires uniform p across elements
  MatrixXd Vq = mesh->Vq;
  MatrixXd xq(Vq.rows(),mesh->K);
  MatrixXd yq(Vq.rows(),mesh->K);
  MatrixXd zq(Vq.rows(),mesh->K);

  // interpolate to cubature points
  for (int e = 0; e < mesh->K; ++e){
    xq.col(e) = Vq*(mesh->x.col(e));
    yq.col(e) = Vq*(mesh->y.col(e));
    zq.col(e) = Vq*(mesh->z.col(e));
  }

  mesh->xq = xq;
  mesh->yq = yq;
  mesh->zq = zq;

  // interp deriv matrices to quadrature
  MatrixXd Vrqtmp,Vsqtmp,Vtqtmp;
  GradVandermonde3D(p_N,mesh->rq,mesh->sq,mesh->tq,Vrqtmp,Vsqtmp,Vtqtmp);
  MatrixXd Vrq = mrdivide(Vrqtmp,mesh->V);
  MatrixXd Vsq = mrdivide(Vsqtmp,mesh->V);
  MatrixXd Vtq = mrdivide(Vtqtmp,mesh->V);

  mesh->Vq = Vq;
  mesh->Vrq = Vrq;
  mesh->Vsq = Vsq;
  mesh->Vtq = Vtq;

  MatrixXd rxq,sxq,txq,ryq,syq,tyq,rzq,szq,tzq,Jq;
  rxq.resize(Vq.rows(),mesh->K);  sxq.resize(Vq.rows(),mesh->K);  txq.resize(Vq.rows(),mesh->K);
  ryq.resize(Vq.rows(),mesh->K);  syq.resize(Vq.rows(),mesh->K);  tyq.resize(Vq.rows(),mesh->K);
  rzq.resize(Vq.rows(),mesh->K);  szq.resize(Vq.rows(),mesh->K);  tzq.resize(Vq.rows(),mesh->K);
  Jq.resize(Vq.rows(),mesh->K);
  for (int e = 0; e < mesh->K; ++e){
    MatrixXd vgeo = vgeofacs3d(mesh->x.col(e),
			       mesh->y.col(e),
			       mesh->z.col(e),
			       Vrq,Vsq,Vtq);
    rxq.col(e) = vgeo.col(0);    sxq.col(e) = vgeo.col(1);    txq.col(e) = vgeo.col(2);
    ryq.col(e) = vgeo.col(3);    syq.col(e) = vgeo.col(4);    tyq.col(e) = vgeo.col(5);
    rzq.col(e) = vgeo.col(6);    szq.col(e) = vgeo.col(7);    tzq.col(e) = vgeo.col(8);
    Jq.col(e)  = vgeo.col(9);
  }
  mesh->Jq = Jq;
  mesh->rxq = rxq;
  mesh->ryq = ryq;
  mesh->rzq = rzq;
  mesh->sxq = sxq;
  mesh->syq = syq;
  mesh->szq = szq;
  mesh->txq = txq;
  mesh->tyq = tyq;
  mesh->tzq = tzq;

  double maxJ = Jq.maxCoeff();
  double minJ = Jq.minCoeff();
  int minJelem = -1;
  for (int e = 0; e < mesh->K; ++e){
    if (Jq.col(e).minCoeff() == minJ){
      minJelem = e;
      break;
    }
  }
  int maxJelem = -1;
  for (int e = 0; e < mesh->K; ++e){
    if (Jq.col(e).maxCoeff() == maxJ){
      maxJelem = e;
      break;
    }
  }
  printf("Minimum Jacobian = %g on elem %d, max Jacobian = %g, on elem %d\n",
	 minJ,minJelem,maxJ,maxJelem);

  // interpolate to face cubature points
  MatrixXd Vfq = mesh->Vfq;

  // interp face nodes to face qnodes matrix
  MatrixXd VfqFace = mesh->VfqFace;
  MatrixXd nxq,nyq,nzq,sJq,Jfq;
  nxq.resize(mesh->rfq.rows(),mesh->K);
  nyq.resize(mesh->rfq.rows(),mesh->K);
  nzq.resize(mesh->rfq.rows(),mesh->K);
  sJq.resize(mesh->rfq.rows(),mesh->K);
  Jfq.resize(mesh->rfq.rows(),mesh->K);
  MatrixXd Vrftmp,Vsftmp,Vtftmp;
  GradVandermonde3D(p_N,mesh->rfq,mesh->sfq,mesh->tfq,Vrftmp,Vsftmp,Vtftmp);
  MatrixXd Vrf = mrdivide(Vrftmp,mesh->V);
  MatrixXd Vsf = mrdivide(Vsftmp,mesh->V);
  MatrixXd Vtf = mrdivide(Vtftmp,mesh->V);

  for (int e = 0; e < mesh->K; ++e){
    MatrixXd sgeo = sgeofacs3d(mesh->x.col(e),
			       mesh->y.col(e),
			       mesh->z.col(e),
			       Vrf,Vsf,Vtf);
    nxq.col(e) = sgeo.col(0);
    nyq.col(e) = sgeo.col(1);
    nzq.col(e) = sgeo.col(2);
    sJq.col(e) = sgeo.col(3);
  }

  mesh->nxq = nxq;
  mesh->nyq = nyq;
  mesh->nzq = nzq;
  mesh->sJq = sJq;

  const int edgenum[6][2] = { {0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3} };
  double hMax = 0.0;
  for (int e = 0; e < mesh->K; ++e){
    double hK1 = Jq.col(e).maxCoeff() / sJq.col(e).maxCoeff();
    double hK2 = 0; // max edge length
    for (int edge = 0; edge < 6; ++edge){
      double dx = mesh->GX(e,edgenum[edge][0])-mesh->GX(e,edgenum[edge][1]);
      double dy = mesh->GY(e,edgenum[edge][0])-mesh->GY(e,edgenum[edge][1]);
      double dz = mesh->GZ(e,edgenum[edge][0])-mesh->GZ(e,edgenum[edge][1]);
      double tmp = sqrt(dx*dx + dy*dy + dz*dz);
      hK2 = max(hK2,tmp);
    }
    double hK3 = pow(Jq.col(e).maxCoeff(),1.0/3.0); // J ~ h^d
    double hK = max(hK1,max(hK2,hK3));
    //double hK = hK2;

    hMax = max(hK,hMax); // take max over entire mesh
  }
  mesh->hMax = hMax;
}
