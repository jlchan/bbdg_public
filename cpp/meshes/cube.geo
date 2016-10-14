
// PARAMETERS

ax = .5;
ay = .5;
az = .5;

Mesh.CharacteristicLengthMax = 10;
Mesh.CharacteristicLengthFactor = .25;
Mesh.Optimize = 1;

BOUNDARY_X = 100;
BOUNDARY_Y = 101;
BOUNDARY_Z = 102;
DOMAIN = 200;


// MESH REFINEMENT AROUND A POINT
// Usefull reference: tutorial 10 of gmsh reference manual

pointAttractor = newp; Point(pointAttractor) = {0, 0, 0};

Field[1] = Attractor;
Field[1].NodesList = {pointAttractor};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = Mesh.CharacteristicLengthMax*10;
Field[2].LcMax = Mesh.CharacteristicLengthMax*10;
Field[2].DistMin = 0.15;
Field[2].DistMax = 0.5;

Background Field = 2;


// GEOMETRY

// Points
ptX = {-ax, ax};
ptY = {-ay, ay};
ptZ = {-az, az};
pIter=0;
For i In {0:1}
  For j In {0:1}
    For k In {0:1}
      p[pIter] = newp; Point(p[pIter]) = {ptX[i], ptY[j], ptZ[k]};
//      p[pIter] = newp; Point(p[pIter]) = {ptX[i], ptY[j], ptZ[k], 1000};
      pIter++;
    EndFor
  EndFor
EndFor

// Lines
lxIter=0;
lyIter=0;
lzIter=0;
For j In {0:1}
  For k In {0:1}
    pNum1 = k + 2*j;
    pNum2 = pNum1+4;
    lx[lxIter] = newl; Line(lx[lxIter]) = {p[pNum1], p[pNum2]};
    lxIter++;
  EndFor
EndFor
For i In {0:1}
  For k In {0:1}
    pNum1 = k + 4*i;
    pNum2 = pNum1+2;
    ly[lyIter] = newl; Line(ly[lyIter]) = {p[pNum1], p[pNum2]};
    lyIter++;
  EndFor
EndFor
For i In {0:1}
  For j In {0:1}
    pNum1 = 2*j + 4*i;
    pNum2 = pNum1+1;
    lz[lzIter] = newl; Line(lz[lzIter]) = {p[pNum1], p[pNum2]};
    lzIter++;
  EndFor
EndFor

// Surfaces
For i In {0:1}
  lyNum1 = 2*i;
  lzNum1 = 2*i;
  lyNum2 = lyNum1 + 1;
  lzNum2 = lzNum1 + 1;
  sxLL = newll; Line Loop(sxLL) = {-lz[lzNum1], ly[lyNum1], lz[lzNum2], -ly[lyNum2]};
  sx[i] = news; Plane Surface(sx[i]) = {sxLL};
EndFor
For j In {0:1}
  lxNum1 = 2*j;
  lzNum1 =   j;
  lxNum2 = lxNum1 + 1;
  lzNum2 = lzNum1 + 2;
  syLL = newll; Line Loop(syLL) = {-lx[lxNum1], lz[lzNum1], lx[lxNum2], -lz[lzNum2]};
  sy[j] = news; Plane Surface(sy[j]) = {syLL};
EndFor
For k In {0:1}
  lxNum1 = k;
  lyNum1 = k;
  lxNum2 = lxNum1 + 2;
  lyNum2 = lyNum1 + 2;
  szLL = newll; Line Loop(szLL) = {-ly[lyNum1], lx[lxNum1], ly[lyNum2], -lx[lxNum2]};
  sz[k] = news; Plane Surface(sz[k]) = {szLL};
EndFor

// Volumes
vSL = newsl; Surface Loop(vSL) = {-sx[0], -sy[0], -sz[0], sx[1], sy[1], sz[1]};
v = news; Volume(v) = {vSL};


// PHYSICAL REGIONS

Physical Surface(BOUNDARY_X) = {-sx[0], sx[1]};
Physical Surface(BOUNDARY_Y) = {-sy[0], sy[1]};
Physical Surface(BOUNDARY_Z) = {-sz[0], sz[1]};
Physical Volume(DOMAIN) = {v};

