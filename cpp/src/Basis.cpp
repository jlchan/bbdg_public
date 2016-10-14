#include <occa.hpp>
#include "Basis.h"
#include "NodeData.h"

// ================== test eigen ================

void test_solve(){
  MatrixXd A(3,3);
  A(0,0) = 1.0;  A(0,1) = 1.0;  A(0,2) = 0.0;
  A(1,0) = 0.0;  A(1,1) = 2.0;  A(1,2) = 2.0;
  A(2,0) = 0.0;  A(2,1) = 0.0;  A(2,2) = 3.0;
  MatrixXd b(3,1);
  b(0,0) = 1.0;
  b(1,0) = 2.0;
  b(2,0) = 3.0;

  MatrixXd x = mldivide(A,b);
  cout << "matrix = " << A << endl;
  cout << "rhs = " << b << endl;
  cout << "solution = " << x << endl;
  MatrixXd err = (A*x - b);
  double relative_error = err.norm() / b.norm(); // norm() is L2 norm
  printf("Err = %f, rel err = %f\n", err.norm(),relative_error);

  printf("mrdivide test:\n");
  MatrixXd bt = b.transpose();
  x = mrdivide(bt,A);
  cout << "matrix = " << A << endl;
  cout << "rhs = " << b << endl;
  cout << "solution = " << x << endl;
  //cout << "Err is " << err << ", while relative error is: " << relative_error << endl;
}

void test_basis(){
  int N = 2;
  VectorXd r(N+1);
  for (int i = 0; i < N+1; ++i){
    r(i) = -1.0 + 2.0*i/N;
  }

  VectorXd P = JacobiP(r, 0.0, 0.0, N);
  VectorXd dP = GradJacobiP(r, 0.0, 0.0, N);
  //  cout << "P = " << P << ", gradP = " << dP << endl;

  MatrixXd V2 = Vandermonde2D(N,r,r);
  cout << "V2 = " << V2 << endl;

  VectorXd dP3 = Simplex3DP(r,r,r,0,1,2);
  //cout << "dP3 = " << dP3 << endl;

  VectorXd pr,ps,pt;
  GradSimplex3DP(r,r,r,1,1,1,pr,ps,pt);
  cout << "pr = " << pr << ", ps = " << ps << ", pt = " << pt << endl;

  VectorXd a,b,c;
  rsttoabc(r,r,r,a,b,c);
  //cout << "a = " << a << ", b = " << b << ", c = " << c << endl;

  MatrixXd V = Vandermonde3D(N,r,r,r);
  MatrixXd Vr,Vs,Vt;
  GradVandermonde3D(N,r,r,r,Vr,Vs,Vt);
  //  cout << "V = " << V << endl;
  //  cout << "Vr = " << Vr << endl;
  //  cout << "Vs = " << Vs << endl;
  //  cout << "Vt = " << Vt << endl;

  MatrixXd VBtri = BernTri(N,r,r);
  MatrixXd VBtet = BernTet(N,r,r,r);
  //cout << "VBtri = " << endl << VBtri << endl;
  //cout << "VBtet = " << endl << VBtet << endl;

  MatrixXd VB1,VB2,VB3,VB4;
  GradBernTet(N,r,r,r,VB1,VB2,VB3,VB4);
  cout << "VB1 = " << endl << VB1 << endl;
  cout << "VB2 = " << endl << VB2 << endl;
  cout << "VB3 = " << endl << VB3 << endl;
  cout << "VB4 = " << endl << VB4 << endl;
}

// =================== begin matlab codes =======================

// d = degree of basis
VectorXd JacobiP(VectorXd x, dfloat alpha, dfloat beta, int d){
  int Nx = x.rows();
  //cout << "Nx = " << Nx << ", x = " << x << endl;
  VectorXd P(Nx,1);
  MatrixXd PL(Nx,d+1);

  // initial P_0 and P_1
  double gamma0 = pow(2,alpha+beta+1)/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);
  for(int i=0; i<Nx; i++)
    PL(i,0) = 1.0/sqrt(gamma0);
  if (d==0)
    return PL;

  double gamma1 = (alpha+1.)*(beta+1.)/(alpha+beta+3.)*gamma0;
  for(int i=0; i<Nx; i++)
    PL(i,1) = ((alpha+beta+2)*x(i)/2.0 + (alpha-beta)/2)/sqrt(gamma1);
  if (d==1){
    for(int i=0; i<Nx; i++)
      P(i) = PL(i,1);
    return P;
  }

  // repeat value in recurrence
  double aold = 2.0/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

  // recurrence
  for(int i=1; i <= d-1; i++){
    double h1 = 2*i + alpha + beta;
    double anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3) );
    double bnew = - (alpha*alpha-beta*beta)/h1/(h1 + 2);
    //cout << "anew = " << anew << ", bnew = " << bnew << endl;
    for(int j=1; j <= Nx; j++){
      PL(j-1,i-1+2) = 1.0/anew*(-aold*PL(j-1,i-1) + (x(j-1) - bnew)*PL(j-1,i-1+1) );
    }
    aold = anew;
  }

  for(int i = 0; i < Nx; i++){
    P(i) = PL(i,d);
  }
  return P;

}

VectorXd GradJacobiP(VectorXd x, double alpha, double beta, int p){

  /// function [dP] = gradJacobiP(z, alpha, beta, N);
  /// Purpose: Evaluate the derivative of the orthonormal Jacobi
  ///   polynomial of type (alpha,beta)>-1, at points x
  ///          for order N and returns dP[1:length(xp))]
  int Nx = x.rows();
  VectorXd dP(Nx);

  if(p == 0)
    dP = 0.0*dP;
  else{
    dP = sqrt(p*(p+alpha+beta+1))*JacobiP(x,alpha+1,beta+1, p-1);
  }
  return dP;
}

// ================= 2D routines ===================

void rstoab(VectorXd r, VectorXd s,
	    VectorXd &a,VectorXd &b){

  int Npts = r.rows();
  a.resize(Npts);
  b.resize(Npts);
  double tol = 1e-8;

  for(int n=0; n<Npts; ++n){
    if(fabs(s(n)-1)>tol)
      a(n) = 2*(1+r(n))/(1-s(n))-1;
    else
      a(n) = -1;
    b(n) = s(n);
  }
}

VectorXd Simplex2DP(VectorXd a,VectorXd b,int i, int j){

  VectorXd h1 = JacobiP(a,0,0,i);
  VectorXd h2 = JacobiP(b,2*i+1,0,j);

  VectorXd bpow = Eigen::pow(1-b.array(),i);

  VectorXd P = sqrt(2)*h1.array()*h2.array()*bpow.array();
  return P;

}

MatrixXd Vandermonde2D(int N, VectorXd r, VectorXd s){
  int Np = (N+1)*(N+2)/2;
  MatrixXd Vout(r.rows(),Np);

  VectorXd a,b,c;
  rstoab(r,s,a,b);
  int sk = 0;
  for (int i = 0; i <= N; ++i){
    for(int j = 0; j <= N - i; ++j){
      Vout.col(sk) = Simplex2DP(a,b,i,j);
      sk = sk+1;
    }
  }
  return Vout;

}

// ========================== 3D ===========================

void rsttoabc(VectorXd r, VectorXd s, VectorXd t,
	      VectorXd &a,VectorXd &b,VectorXd &c){

  int Np = r.rows();
  a.resize(Np);
  b.resize(Np);
  c.resize(Np);
  double tol = 1e-8;
  for(int n=0; n<Np; ++n){
    if(fabs(s(n)+t(n))>tol)
      a(n) = 2*(1+r(n))/(-s(n)-t(n))-1;
    else
      a(n) = -1;

    if(fabs(t(n)-1)>tol)
      b(n) = 2*(1+s(n))/(1-t(n))-1;
    else
      b(n) = -1;

    c(n) = t(n);
  }
}

VectorXd Simplex3DP(VectorXd a, VectorXd b, VectorXd c,
		    int i, int j, int k){

  // function [P] = Simplex3DP(a,b,c,i,j,k);
  // Purpose : Evaluate 3D orthonormal polynomial
  //           on simplex at (a,b,c) of order (i,j,k).

  VectorXd h1 = JacobiP(a,0,0,i);
  VectorXd h2 = JacobiP(b,2*i+1,0,j);
  VectorXd h3 = JacobiP(c,2*(i+j)+2,0,k);

  VectorXd bpow = Eigen::pow(1-b.array(),i);
  VectorXd cpow = Eigen::pow(1-c.array(),i+j);

  VectorXd P = 2*sqrt(2)*h1.array()*h2.array()*bpow.array()*h3.array()*cpow.array();
  return P;
}

void GradSimplex3DP(VectorXd a, VectorXd b, VectorXd c,
		    int id, int jd, int kd,
		    VectorXd &V3Dr, VectorXd &V3Ds, VectorXd &V3Dt){

  // function [V3Dr, V3Ds, V3Dt] = GradSimplex3DP(a,b,c,id,jd,kd)
  // Purpose: Return the derivatives of the modal basis (id,jd,kd)
  //          on the 3D simplex at (a,b,c)

  VectorXd fa,dfa,gb,dgb,hc,dhc;
  fa = JacobiP(a,0,0,id);           dfa = GradJacobiP(a,0,0,id);
  gb = JacobiP(b,2*id+1,0,jd);      dgb = GradJacobiP(b,2*id+1,0,jd);
  hc = JacobiP(c,2*(id+jd)+2,0,kd); dhc = GradJacobiP(c,2*(id+jd)+2,0,kd);

  // r-derivative
  int Np = a.rows();
  V3Dr.resize(Np);
  V3Dr.array() = dfa.array()*(gb.array()*hc.array());
  if(id>0)
    V3Dr.array() *= Eigen::pow(0.5*(1-b.array()),id-1);

  if(id+jd>0)
    V3Dr.array() *= Eigen::pow(0.5*(1-c.array()),id+jd-1);

  // s-derivative
  V3Ds.resize(Np);
  V3Ds.array() = 0.5*(1+a.array())*V3Dr.array();
  VectorXd tmp = dgb.array()*Eigen::pow(0.5*(1-b.array()),id);
  if(id > 0)
    tmp.array() += (-0.5*id)*(gb.array()*Eigen::pow(0.5*(1-b.array()),id-1));
  if(id+jd>0)
    tmp.array() *= Eigen::pow(0.5*(1-c.array()),id+jd-1);
  tmp.array() = fa.array()*(tmp.array()*hc.array());
  V3Ds.array() += tmp.array();

  // t-derivative
  V3Dt.resize(Np);
  V3Dt.array() = 0.5*(1+a.array())*V3Dr.array() + 0.5*(1+b.array())*tmp.array();
  tmp.array() = dhc.array()*Eigen::pow(0.5*(1-c.array()),id+jd);
  if(id+jd>0)
    tmp.array() += -0.5*(id+jd)*(hc.array()*Eigen::pow(0.5*(1-c.array()),id+jd-1));
  tmp.array() = fa.array()*(gb.array()*tmp.array());
  tmp.array() *= Eigen::pow(0.5*(1-b.array()),id);
  V3Dt.array() = V3Dt.array()+tmp.array();

  // normalize
  V3Dr = V3Dr*pow(2,2*id+jd+1.5);
  V3Ds = V3Ds*pow(2,2*id+jd+1.5);
  V3Dt = V3Dt*pow(2,2*id+jd+1.5);

}

MatrixXd Vandermonde3D(int N, VectorXd r, VectorXd s, VectorXd t){
  /// function [Vout] = Voutandermonde(N,r,s,t,Vout)
  /// Purpose : Initialize the gradient of the modal basis (i,j,k)
  ///at (r,s,t) at order p

  int Np = (N+1)*(N+2)*(N+3)/6;
  MatrixXd Vout(r.rows(),Np);

  /// find tensor-product coordinates
  VectorXd a, b, c;
  rsttoabc(r,s,t,a,b,c);

  /// initialize matrices
  int sk=0;
  for(int i=0; i<=N; ++i){
    for(int j=0; j<=N-i; ++j){
      for(int k=0; k<=N-i-j; ++k){
	Vout.col(sk) = Simplex3DP(a,b,c,i,j,k);
        sk++;
      }
    }
  }
  return Vout;
}


void GradVandermonde3D(int N, VectorXd r, VectorXd s, VectorXd t,
		       MatrixXd &V3Dr, MatrixXd &V3Ds, MatrixXd &V3Dt){

  int Npts = r.rows();
  int Np = (N+1)*(N+2)*(N+3)/6;
  V3Dr.resize(Npts,Np);
  V3Ds.resize(Npts,Np);
  V3Dt.resize(Npts,Np);

  // find tensor-product coordinates
  VectorXd a,b,c;
  rsttoabc(r,s,t,a,b,c);

  // Initialize matrices
  int sk = 0;
  for (int i = 0; i <= N; ++i){
    for (int j = 0; j <= N-i; ++j){
      for (int k = 0; k <= N-i-j; ++k){
	VectorXd Vr,Vs,Vt;
	GradSimplex3DP(a,b,c,i,j,k,Vr,Vs,Vt);
	V3Dr.col(sk) = Vr;
	V3Ds.col(sk) = Vs;
	V3Dt.col(sk) = Vt;
	sk = sk+1;
      }
    }
  }
}

// =============== Hierarchical basis ==================

void barytors(VectorXd L1, VectorXd L2, VectorXd L3, VectorXd &r,VectorXd &s){
  VectorXd v1(2); v1 << -1,-1;
  VectorXd v2(2); v2 << 1,-1;
  VectorXd v3(2); v3 << -1,1;

  int Npts = L1.rows();
  r.resize(Npts);
  s.resize(Npts);
  for (int i = 0; i < Npts; ++i){
    VectorXd XY = v1*L1(i) + v2*L2(i) + v3*L3(i);
    r(i) = XY(1);
    s(i) = XY(2);
  }
}

void VandermondeHier(int N, VectorXd r, VectorXd s, VectorXd t,
		     MatrixXd &vertV, MatrixXd &edgeV, MatrixXd &faceV){

  int Np = (N+1)*(N+2)*(N+3)/6;
  int Npts = r.rows();
  vertV.resize(Npts,4);
  edgeV.resize(Npts,6*(N-1));
  faceV.resize(Npts,4*(N-2)*(N-1)/2);

  // vertex functions
  vertV.col(0) = -(1.0 + r.array() + s.array() + t.array())/2.0;
  vertV.col(1) = (1.0 + r.array())/2.0;
  vertV.col(2) = (1.0 + s.array())/2.0;
  vertV.col(3) = (1.0 + t.array())/2.0;

  MatrixXi edges(6,2);
  edges << 0,3,
    1, 3,
    2, 3,
    0, 1,
    1, 2,
    2, 0;

  // edge basis functions
  int sk = 0;
  for (int edge =0; edge < 6; ++edge){
    for (int i = 0; i <= N-2; ++i){
      int vertex1 = edges(edge,0);
      int vertex2 = edges(edge,1);
      VectorXd V1,V2;
      V1 = vertV.col(vertex1);
      V2 = vertV.col(vertex2);
      VectorXd JP = JacobiP(V1.array()-V2.array(),0,0,i);
      edgeV.col(sk) = V1.array()*V2.array()*JacobiP(V1.array()-V2.array(),0,0,i).array();
      ++sk;
    }
  }

  MatrixXi faces(4,3);
  faces << 0,1,3,
    1,2,3,
    2,0,3,
    0,1,2;

  // face basis functions
  sk = 0;
  VectorXd rr,ss;
  for (int f = 0; f < 4; ++f){

    int vertex1 = faces(f,0);
    int vertex2 = faces(f,1);
    int vertex3 = faces(f,2);

    VectorXd V1 = vertV.col(vertex1);
    VectorXd V2 = vertV.col(vertex2);
    VectorXd V3 = vertV.col(vertex3);

    barytors(V1,V2,V3,rr,ss); // convert affine coords to r,s
    MatrixXd Vface = Vandermonde2D(N-3,rr,ss);
    for (int i = 0; i < Vface.cols(); ++i){
      faceV.col(sk) = V1.array()*V2.array()*V3.array()*Vface.col(i).array();
      ++sk;
    }
  }
}

// hierarchical basis over all vertex/edges/faces - can be used to determine nodes
MatrixXd VandermondeGHsurf(int N, int Npsurf, VectorXd r, VectorXd s, VectorXd t){
  int Np = (N+1)*(N+2)*(N+3)/6;
  int Npts = r.rows();
  MatrixXd V(Npts,Npsurf);

  // vertex functions
  V.col(0) = -(1.0 + r.array() + s.array() + t.array())/2.0;
  V.col(1) = (1.0 + r.array())/2.0;
  V.col(2) = (1.0 + s.array())/2.0;
  V.col(3) = (1.0 + t.array())/2.0;

  MatrixXi edges(6,2);
  edges << 0,3,
    1, 3,
    2, 3,
    0, 1,
    1, 2,
    2, 0;

  // edge basis functions
  int sk = 4;
  for (int edge =0; edge < 6; ++edge){
    for (int i = 0; i <= N-2; ++i){
      int vertex1 = edges(edge,0);
      int vertex2 = edges(edge,1);
      VectorXd V1,V2;
      V1 = V.col(vertex1);
      V2 = V.col(vertex2);
      VectorXd JP = JacobiP(V1.array()-V2.array(),0,0,i);
      V.col(sk) = V1.array()*V2.array()*JacobiP(V1.array()-V2.array(),0,0,i).array();
      ++sk;
    }
  }

  MatrixXi faces(4,3);
  faces << 0,1,3,
    1,2,3,
    2,0,3,
    0,1,2;

  // face basis functions
  VectorXd rr,ss;
  for (int f = 0; f < 4; ++f){

    int vertex1 = faces(f,0);
    int vertex2 = faces(f,1);
    int vertex3 = faces(f,2);

    VectorXd V1 = V.col(vertex1);
    VectorXd V2 = V.col(vertex2);
    VectorXd V3 = V.col(vertex3);

    barytors(V1,V2,V3,rr,ss); // convert affine coords to r,s
    MatrixXd Vface = Vandermonde2D(N-3,rr,ss);
    for (int i = 0; i < Vface.cols(); ++i){
      V.col(sk) = V1.array()*V2.array()*V3.array()*Vface.col(i).array();
      ++sk;
    }
  }
  return V;
}

// =============== Bernstein basis ================

unsigned int factorial_ratio(int n1,int n2){
  unsigned int val = 1;
  for (int i = n1; i > n2; --i){
    val *= i;
  }
  return val;
}

unsigned int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

unsigned int nchoosek(int n,int k){
  return factorial(n) / (factorial(n-k) * factorial(k));
}

MatrixXd Bern1D(int N, VectorXd r){
  VectorXd x = .5*(1.0+r.array());
  int Np = (N+1);
  MatrixXd V(r.rows(),Np); V.fill(0.0);
  int sk = 0;
  for (int i = 0; i < Np; ++i){
    V.col(i) = nchoosek(N,i) * x.array().pow(i) * (1.0 - x.array()).array().pow(N-i);
  }
  return V;

}

MatrixXd BernTri(int N, VectorXd r, VectorXd s){
  VectorXd L1 = -(r.array()+s.array())/2;
  VectorXd L2 = (1+r.array())/2;
  VectorXd L3 = (1+s.array())/2;

  int Np = (N+1)*(N+2)/2;
  MatrixXd V(r.rows(),Np);
  int sk = 0;
  for (int k = 0; k <= N; ++k){
    for (int j = 0; j <= N-k; ++j){
      int i = N-j-k;
      double C = (double) factorial(N)/(factorial(i)*factorial(j)*factorial(k));
      VectorXd L1i = Eigen::pow(L1.array(),i);
      VectorXd L2j = Eigen::pow(L2.array(),j);
      VectorXd L3k = Eigen::pow(L3.array(),k);
      V.col(sk) = C*L1i.array()*L2j.array()*L3k.array();
      sk = sk + 1;
    }
  }
  return V;
}

MatrixXd BernTet(int N, VectorXd r, VectorXd s, VectorXd t){

  int Np = (N+1)*(N+2)*(N+3)/6;
  MatrixXd V(r.rows(),Np);

  VectorXd L1 = -(1+r.array()+s.array()+t.array())/2;
  VectorXd L2 = (1+r.array())/2;
  VectorXd L3 = (1+s.array())/2;
  VectorXd L4 = (1+t.array())/2;

  int sk = 0;
  for (int l = 0; l <= N; ++l){
    for (int k = 0; k <= N-l; ++k){
      for (int j = 0; j <= N-l-k; ++j){
	int i = N-j-k-l;
	double C = (double) factorial(N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
	VectorXd L1i = Eigen::pow(L1.array(),i);
	VectorXd L2j = Eigen::pow(L2.array(),j);
	VectorXd L3k = Eigen::pow(L3.array(),k);
	VectorXd L4l = Eigen::pow(L4.array(),l);
	V.col(sk) = C*L1i.array()*L2j.array()*L3k.array()*L4l.array();
	++sk;
      }
    }
  }

  return V;
}

void GradBernTet(int N, VectorXd r, VectorXd s, VectorXd t,
		 MatrixXd &V1, MatrixXd &V2, MatrixXd &V3, MatrixXd &V4){

  int Np = (N+1)*(N+2)*(N+3)/6;

  VectorXd L1 = -(1+r.array()+s.array()+t.array())/2;
  VectorXd L2 = (1+r.array())/2;
  VectorXd L3 = (1+s.array())/2;
  VectorXd L4 = (1+t.array())/2;

  V1.resize(r.rows(),Np);
  V2.resize(r.rows(),Np);
  V3.resize(r.rows(),Np);
  V4.resize(r.rows(),Np);
  int sk = 0;
  for (int l = 0; l <= N; ++l){
    for (int k = 0; k <= N-l; ++k){
      for (int j = 0; j <= N-l-k; ++j){
	int i = N-j-k-l;
	double C = (double) factorial(N)/(factorial(i)*factorial(j)*factorial(k)*factorial(l));
	VectorXd L1i = Eigen::pow(L1.array(),i);
	VectorXd L2j = Eigen::pow(L2.array(),j);
	VectorXd L3k = Eigen::pow(L3.array(),k);
	VectorXd L4l = Eigen::pow(L4.array(),l);

	VectorXd L1im = Eigen::pow(L1.array(),i-1);
	VectorXd L2jm = Eigen::pow(L2.array(),j-1);
	VectorXd L3km = Eigen::pow(L3.array(),k-1);
	VectorXd L4lm = Eigen::pow(L4.array(),l-1);

	VectorXd dL1 = C*i*L1im.array()*L2j.array()*L3k.array()*L4l.array();
	VectorXd dL2 = C*j*L1i.array()*L2jm.array()*L3k.array()*L4l.array();
	VectorXd dL3 = C*k*L1i.array()*L2j.array()*L3km.array()*L4l.array();
	VectorXd dL4 = C*l*L1i.array()*L2j.array()*L3k.array()*L4lm.array();
	if (i==0)
	  dL1.fill(0.0);
	if (j==0)
	  dL2.fill(0.0);
	if (k ==0)
	  dL3.fill(0.0);
	if (l ==0)
	  dL4.fill(0.0);

	V1.col(sk) = dL1;
	V2.col(sk) = dL2;
	V3.col(sk) = dL3;
	V4.col(sk) = dL4;

	++sk;
      }
    }
  }
}

// ===================== geofacs routines =======================

// (v/s)geofacs = array of volume, surface geofacs
MatrixXd vgeofacs3d(VectorXd x,VectorXd y,VectorXd z,
		    MatrixXd Dr,MatrixXd Ds,MatrixXd Dt){

  int Npts = Dr.rows();
  MatrixXd vgeofacs(Npts,10); // rx,sx,tx,ry,sy,ty,rz,sz,tz,J

  VectorXd xr = Dr*x;  VectorXd xs = Ds*x;  VectorXd xt = Dt*x;
  VectorXd yr = Dr*y;  VectorXd ys = Ds*y;  VectorXd yt = Dt*y;
  VectorXd zr = Dr*z;  VectorXd zs = Ds*z;  VectorXd zt = Dt*z;

  VectorXd J =
    xr.array()*(ys.array()*zt.array()-zs.array()*yt.array()) -
    yr.array()*(xs.array()*zt.array()-zs.array()*xt.array()) +
    zr.array()*(xs.array()*yt.array()-ys.array()*xt.array());

  VectorXd rx =  (ys.array()*zt.array() - zs.array()*yt.array())/J.array();
  VectorXd ry = -(xs.array()*zt.array() - zs.array()*xt.array())/J.array();
  VectorXd rz =  (xs.array()*yt.array() - ys.array()*xt.array())/J.array();
  VectorXd sx = -(yr.array()*zt.array() - zr.array()*yt.array())/J.array();
  VectorXd sy =  (xr.array()*zt.array() - zr.array()*xt.array())/J.array();
  VectorXd sz = -(xr.array()*yt.array() - yr.array()*xt.array())/J.array();
  VectorXd tx =  (yr.array()*zs.array() - zr.array()*ys.array())/J.array();
  VectorXd ty = -(xr.array()*zs.array() - zr.array()*xs.array())/J.array();
  VectorXd tz =  (xr.array()*ys.array() - yr.array()*xs.array())/J.array();

  vgeofacs.col(0) = rx;
  vgeofacs.col(1) = sx;
  vgeofacs.col(2) = tx;
  vgeofacs.col(3) = ry;
  vgeofacs.col(4) = sy;
  vgeofacs.col(5) = ty;
  vgeofacs.col(6) = rz;
  vgeofacs.col(7) = sz;
  vgeofacs.col(8) = tz;
  vgeofacs.col(9) = J;
  return vgeofacs;
}

// D(rst)f = matrix mapping nodal values to derivatives at face cubature points
MatrixXd sgeofacs3d(VectorXd x,VectorXd y,VectorXd z,
		    MatrixXd Drf,MatrixXd Dsf,MatrixXd Dtf){
  MatrixXd vgeo = vgeofacs3d(x,y,z,Drf,Dsf,Dtf);
  VectorXd rx = vgeo.col(0);
  VectorXd sx = vgeo.col(1);
  VectorXd tx = vgeo.col(2);
  VectorXd ry = vgeo.col(3);
  VectorXd sy = vgeo.col(4);
  VectorXd ty = vgeo.col(5);
  VectorXd rz = vgeo.col(6);
  VectorXd sz = vgeo.col(7);
  VectorXd tz = vgeo.col(8);
  VectorXd J  = vgeo.col(9);

  int Npts = Drf.rows();
  MatrixXd sgeofacs(Npts,4); // nx,ny,nz,sJ

  int Nfpts = Npts/4; // assume each face has the same quadrature
  VectorXd nx(Npts);
  VectorXd ny(Npts);
  VectorXd nz(Npts);
  for (int i = 0; i < Nfpts; ++i){
    // face 1
    int id = i;
    nx(id) = -tx(id);
    ny(id) = -ty(id);
    nz(id) = -tz(id);

    // face 2
    id += Nfpts;
    nx(id) = -sx(id);
    ny(id) = -sy(id);
    nz(id) = -sz(id);

    // face 3
    id += Nfpts;
    nx(id) = rx(id) + sx(id) + tx(id);
    ny(id) = ry(id) + sy(id) + ty(id);
    nz(id) = rz(id) + sz(id) + tz(id);

    // face 4
    id += Nfpts;
    nx(id) = -rx(id);
    ny(id) = -ry(id);
    nz(id) = -rz(id);
  }
  VectorXd sJ =
    (nx.array()*nx.array() +
     ny.array()*ny.array() +
     nz.array()*nz.array()).array().sqrt();

  // normalize and scale sJ
  nx.array()/=sJ.array();
  ny.array()/=sJ.array();
  nz.array()/=sJ.array();
  sJ.array() *= J.array();

  sgeofacs.col(0) = nx;
  sgeofacs.col(1) = ny;
  sgeofacs.col(2) = nz;
  sgeofacs.col(3) = sJ;
  return sgeofacs;

}


// ===================== tabulated nodes ========================

/*
void tet_cubature_duffy(int N, VectorXd &a, VectorXd &wa,
			VectorXd &b,  VectorXd &wb,
			VectorXd &c,  VectorXd &wc){

  int Np = (N+1);
  a.resize(Np); wa.resize(Np);
  b.resize(Np); wb.resize(Np);
  c.resize(Np); wc.resize(Np);

#define loadnodes_duffy(M)						\
  a(i) = p_r1D_N##M[i]; wa(i) = p_s_N##M[i]; t(i) = p_t_N##M[i];

}
*/

// requires alpha, beta = integers
void JacobiGQ(int N, int alpha_int, int beta_int, VectorXd &r, VectorXd &w){

  double alpha = (double) alpha_int;
  double beta = (double) beta_int;

  int Np = (N+1);
  r.resize(Np);
  w.resize(Np);

  if (N==0){
    r(0) = -(alpha + beta)/(alpha + beta + 2.0);
    w(0) = 2.0;
  }
  VectorXd h1(Np);
  for (int i = 0; i < Np; ++i){
    h1(i) = 2*i+alpha+beta;
  }
  MatrixXd J(Np,Np); J.fill(0.0);
  for (int i = 0; i < Np; ++i){
    J(i,i) = -.5*(alpha*alpha-beta*beta) / (h1(i)+2.0) / h1(i);
    if (i < N){
      J(i,i+1) = 2.0 / (h1(i)+2.0) *
	sqrt((i+1.0)*(i+1.0 + alpha+beta)*
	     (i+1.0 + alpha)* (i+1.0 + beta) /
	     (h1(i)+1.0)/ (h1(i)+3.0));
    }
  }
  double tol = 1e-14;
  if (alpha + beta < tol){
    J(0,0) = 0.0;
  }
  MatrixXd Js = J + J.transpose();

  SelfAdjointEigenSolver<MatrixXd> eig(Js);
  MatrixXd V = eig.eigenvectors();
  r = eig.eigenvalues();

  double gamma_alpha = factorial(alpha-1);
  double gamma_beta = factorial(alpha-1);

  // assumes alpha,beta = int so gamma(x) = factorial(x-1)
  w = V.row(0).array().square() * pow(2.0,alpha+beta+1.0) / (alpha+beta+1) *
    factorial(alpha) * factorial(beta) / factorial(alpha+beta);

  //  cout << "r = " << r << endl;
  //  cout << "w = " << w << endl;
}



void Nodes3D(int N, VectorXd &r, VectorXd &s, VectorXd &t){
  int Np = (N+1)*(N+2)*(N+3)/6;
  r.resize(Np);
  s.resize(Np);
  t.resize(Np);
#define loadnodes(M)				\
  r(i) = p_r_N##M[i]; s(i) = p_s_N##M[i]; t(i) = p_t_N##M[i];

  for (int i = 0; i < Np; ++i){
    switch(N){
    case 1: loadnodes(1); break;
    case 2: loadnodes(2); break;
    case 3: loadnodes(3); break;
    case 4: loadnodes(4); break;
    case 5: loadnodes(5); break;
    case 6: loadnodes(6); break;
    case 7: loadnodes(7); break;
    case 8: loadnodes(8); break;
    case 9: loadnodes(9); break;
    case 10: loadnodes(10); break;
    case 11: loadnodes(11); break;
    case 12: loadnodes(12); break;
    case 13: loadnodes(13); break;
    case 14: loadnodes(14); break;
    case 15: loadnodes(15); break;
    case 16: loadnodes(16); break;
    case 17: loadnodes(17); break;
    case 18: loadnodes(18); break;
    case 19: loadnodes(10); break;
    case 20: loadnodes(20); break;
    }
  }
}

void tet_cubature(int N, VectorXd &rq, VectorXd &sq, VectorXd &tq, VectorXd &wq){

#define loadNq(M)  Nq = p_Nq_N##M;
#define loadq(M)				\
  rq(i) = p_rq_N##M[i]; sq(i) = p_sq_N##M[i]; tq(i) = p_tq_N##M[i]; wq(i) = p_wq_N##M[i];

  int Nq;
  switch(N){
  case 1: loadNq(1); break;
  case 2: loadNq(2); break;
  case 3: loadNq(3); break;
  case 4: loadNq(4); break;
  case 5: loadNq(5); break;
  case 6: loadNq(6); break;
  case 7: loadNq(7); break;
  case 8: loadNq(8); break;
  case 9: loadNq(9); break;
  case 10: loadNq(10); break;
  case 11: loadNq(11); break;
  case 12: loadNq(12); break;
  case 13: loadNq(13); break;
  case 14: loadNq(14); break;
  case 15: loadNq(15); break;
  case 16: loadNq(16); break;
  case 17: loadNq(17); break;
  case 18: loadNq(18); break;
  case 19: loadNq(19); break;
  case 20: loadNq(20); break;
  case 21: loadNq(21); break;
  }

  rq.resize(Nq);  sq.resize(Nq);  tq.resize(Nq);  wq.resize(Nq);
  for (int i = 0; i < Nq; ++i){
    switch(N){
    case 1: loadq(1); break;
    case 2: loadq(2); break;
    case 3: loadq(3); break;
    case 4: loadq(4); break;
    case 5: loadq(5); break;
    case 6: loadq(6); break;
    case 7: loadq(7); break;
    case 8: loadq(8); break;
    case 9: loadq(9); break;
    case 10: loadq(10); break;
    case 11: loadq(11); break;
    case 12: loadq(12); break;
    case 13: loadq(13); break;
    case 14: loadq(14); break;
    case 15: loadq(15); break;
    case 16: loadq(16); break;
    case 17: loadq(17); break;
    case 18: loadq(18); break;
    case 19: loadq(19); break;
    case 20: loadq(20); break;
    case 21: loadq(21); break;
    }
  }
}

void tri_cubature(int N, VectorXd &rfq, VectorXd &sfq, VectorXd &wfq){

#define loadNfq(M)  Nfq = p_Nfq_N##M;
#define loadfq(M)				\
  rfq(i) = p_rfq_N##M[i]; sfq(i) = p_sfq_N##M[i]; wfq(i) = p_wfq_N##M[i];

  int Nfq;
  switch(N){
  case 1: loadNfq(1); break;
  case 2: loadNfq(2); break;
  case 3: loadNfq(3); break;
  case 4: loadNfq(4); break;
  case 5: loadNfq(5); break;
  case 6: loadNfq(6); break;
  case 7: loadNfq(7); break;
  case 8: loadNfq(8); break;
  case 9: loadNfq(9); break;
  case 10: loadNfq(10); break;
  case 11: loadNfq(11); break;
  case 12: loadNfq(12); break;
  case 13: loadNfq(13); break;
  case 14: loadNfq(14); break;
  case 15: loadNfq(15); break;
  case 16: loadNfq(16); break;
  case 17: loadNfq(17); break;
  case 18: loadNfq(18); break;
  case 19: loadNfq(19); break;
  case 20: loadNfq(20); break;
  case 21: loadNfq(21); break;
  }

  rfq.resize(Nfq);  sfq.resize(Nfq);  wfq.resize(Nfq);
  for (int i = 0; i < Nfq; ++i){
    switch(N){
    case 1: loadfq(1); break;
    case 2: loadfq(2); break;
    case 3: loadfq(3); break;
    case 4: loadfq(4); break;
    case 5: loadfq(5); break;
    case 6: loadfq(6); break;
    case 7: loadfq(7); break;
    case 8: loadfq(8); break;
    case 9: loadfq(9); break;
    case 10: loadfq(10); break;
    case 11: loadfq(11); break;
    case 12: loadfq(12); break;
    case 13: loadfq(13); break;
    case 14: loadfq(14); break;
    case 15: loadfq(15); break;
    case 16: loadfq(16); break;
    case 17: loadfq(17); break;
    case 18: loadfq(18); break;
    case 19: loadfq(19); break;
    case 20: loadfq(20); break;
    case 21: loadfq(21); break;
    }
  }
}


// ===================== lin algo subroutines =====================


// backslash
MatrixXd mldivide(MatrixXd &A, MatrixXd &B){
  //return A.colPivHouseholderQr().solve(B);
  return A.fullPivHouseholderQr().solve(B);
}

// A/B
MatrixXd mrdivide(MatrixXd &A, MatrixXd &B){
  MatrixXd At = A.transpose();
  MatrixXd Bt = B.transpose();
  //return (Bt.colPivHouseholderQr().solve(At)).transpose();
  return (Bt.fullPivHouseholderQr().solve(At)).transpose();
}

// extract sub-vector routine - currently unused
VectorXd extract(const VectorXd& full, const VectorXi &ind){
  int num_indices = ind.rows();
  VectorXd target(num_indices);
  for (int i = 0; i < num_indices; i++){
    target(i) = full(ind(i));
  }
  return target;
}

// fixed-bandwidth sparse ids
void get_sparse_ids(MatrixXd A, MatrixXi &cols, MatrixXd &vals){
  int nrows = A.rows();
  int ncols = A.cols();

  double maxVal = A.array().abs().maxCoeff();
  double tol = 1e-5;
  MatrixXi boolA = ((A.array().abs())>tol*maxVal).cast<int>();
  int max_col_vals = (boolA.rowwise().sum()).maxCoeff();
  //printf("max col vals = %d\n",max_col_vals);
  cols.resize(nrows,max_col_vals);
  vals.resize(nrows,max_col_vals);

  cols.fill(0);
  vals.fill(0.0);
  //  cout << boolA << endl;
  //  cout << "max col vals = " << max_col_vals << endl;
  for (int i = 0; i < nrows; ++i){
    int sk = 0;
    for (int j = 0; j < ncols; ++j){
      if (fabs(A(i,j))>tol*maxVal){
	cols(i,sk) = j;
	vals(i,sk) = A(i,j);
	++sk;
      }
    }
  }
  //  cout << "vals = " << endl << vals << endl;
}
