
function buildsource(in_N)

if nargin==0
    in_N=1;
end

Globals2D;

NODETOL = 1e-8;

N = in_N;

Np = (N+1)*(N+2)/2; Nfp = N+1; Nfaces = 3;

[x,y] = Nodes2D(N);
[r,s] = xytors(x,y);

% find all the nodes that lie on each edge
fmask1   = find( abs(s+1) < NODETOL)';
fmask2   = find( abs(r+s) < NODETOL)';
fmask3   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3]';

% vandermonde matrix
V = Vandermonde2D(N, r, s);

% derivative matrices
[Dr, Ds] = Dmatrices2D(N, r, s, V);

%% lift matrix
Emat = zeros(Np, Nfaces*Nfp);

for face=1:Nfaces
  % process face
  if(face==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end;
  if(face==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end;
  if(face==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end;
  if(face==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end;
  
  VFace = Vandermonde2D(N, faceR, faceS);
  massFace = inv(VFace*VFace');

  idr = Fmask(:,face); idc = (face-1)*Nfp+1:face*Nfp;

  Emat(idr, idc) = Emat(idr, idc)+ massFace;
end

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);
keyboard

% DnLift


% correct to 0-index
Fmask = Fmask-1;


fid = fopen(sprintf('dataN%02d.h', N), 'w');

fprintf(fid, sprintf('#ifndef DATAN%02d \n\n', N));
fprintf(fid, sprintf('#define DATAN%02d 1\n\n', N));
fprintf(fid, 'double p_r[%d] = { ', Np);
fprintf(fid, '%17.15g ', r(1));
for n=2:Np
    fprintf(fid, ', %17.15g ', r(n));
end
fprintf(fid, '};\n')

fprintf(fid, 'double p_s[%d] = {', Np);
fprintf(fid, '%17.15g ', s(1));
for n=2:Np
    fprintf(fid, ', %17.15g ', s(n));
end
fprintf(fid, '};\n')

fprintf(fid, 'double p_Dr[%d][%d] = {', Np, Np);
for m=1:Np
    fprintf(fid, '{%17.15g ', Dr(m,1));
    for n=2:Np
        fprintf(fid, ', %17.15g ', Dr(m,n));
    end
    if(m<Np)
        fprintf(fid, '},\n')
    else
        fprintf(fid, '}};\n')
    end
end

fprintf(fid, 'double p_Ds[%d][%d] = {', Np, Np);
for m=1:Np
    fprintf(fid, '{%17.15g ', Ds(m,1));
    for n=2:Np
        fprintf(fid, ', %17.15g ', Ds(m,n));
    end
    if(m<Np)
        fprintf(fid, '},\n')
    else
        fprintf(fid, '}};\n')
    end
end

fprintf(fid, 'double p_LIFT[%d][%d] = {', Np, Nfp*Nfaces);
for m=1:Np
    fprintf(fid, '{%17.15g ', LIFT(m,1));
    for n=2:Nfp*Nfaces
        fprintf(fid, ', %17.15g ', LIFT(m,n));
    end
    if(m<Np)
        fprintf(fid, '},\n')
    else
        fprintf(fid, '}};\n')
    end
end


fprintf(fid, 'int p_Fmask[%d][%d] = {', Nfaces, Nfp);
for m=1:Nfaces
    fprintf(fid, '{%d ', Fmask(1,m));
    for n=2:Nfp
        fprintf(fid, ', %d ', Fmask(n,m));
    end
    if(m<Nfaces)
        fprintf(fid, '},\n')
    else
        fprintf(fid, '}};\n')
    end
end

fprintf(fid, sprintf('\n#endif\n'));
fclose(fid)
