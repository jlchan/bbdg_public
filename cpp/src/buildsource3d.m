% writes header files for nodal bases 

function buildsource3d

for N = 1:7
    write_headers(N);
end

function write_headers(in_N)

if nargin==0
    in_N=1;
end

Globals3D;

NODETOL = 1e-8;

N = in_N;

Np = (N+1)*(N+2)*(N+3)/6; Nfp = (N+1)*(N+2)/2; Nfaces = 4;

[x,y,z] = Nodes3D(N);
[r,s,t] = xyztorst(x,y,z);

% find all the nodes that lie on each edge
fmask1   = find( abs(t+1) < NODETOL)';
fmask2   = find( abs(s+1) < NODETOL)';
fmask3   = find( abs(r+s+t+1) < NODETOL)';
fmask4   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

% vandermonde matrix
V = Vandermonde3D(N, r, s, t);

% derivative matrices
[Dr, Ds, Dt] = Dmatrices3D(N, r, s, t, V);

% interp to quadrature points
[rq sq tq w] = tet_cubature(2*N+1);
Vq = Vandermonde3D(N, rq, sq, tq)/V;

%% lift matrix
Emat = zeros(Np, Nfaces*Nfp);

for face=1:Nfaces
  % process face
  if(face==1); faceR = r(Fmask(:,1)); faceS = s(Fmask(:,1)); end;
  if(face==2); faceR = r(Fmask(:,2)); faceS = t(Fmask(:,2)); end;
  if(face==3); faceR = s(Fmask(:,3)); faceS = t(Fmask(:,3)); end;
  if(face==4); faceR = s(Fmask(:,4)); faceS = t(Fmask(:,4)); end;
  
  VFace = Vandermonde2D(N, faceR, faceS);
  massFace = inv(VFace*VFace');

  idr = Fmask(:,face); 
  idc = (face-1)*Nfp+1:face*Nfp;

  Emat(idr, idc) = Emat(idr, idc)+ massFace;
end

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
LIFT = V*(V'*Emat);

% correct to 0-index
Fmask = Fmask-1;

%%

fid = fopen(sprintf('../include/data3dN0%d.h', N), 'w');

%fprintf(fid, sprintf('#ifndef DATAN%03d \n\n', N));
%fprintf(fid, sprintf('#define DATAN%03d 1\n\n', N));
fprintf(fid, sprintf('#ifndef DATA3dN0%d \n\n', N));
fprintf(fid, sprintf('#define DATA3dN0%d 1\n\n', N));

fprintf(fid, 'dfloat p_r[%d] = { ', Np);
fprintf(fid, '%17.15g ', r(1));
for n=2:Np
    fprintf(fid, ', %17.15g ', r(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_s[%d] = {', Np);
fprintf(fid, '%17.15g ', s(1));
for n=2:Np
    fprintf(fid, ', %17.15g ', s(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_t[%d] = {', Np);
fprintf(fid, '%17.15g ', t(1));
for n=2:Np
    fprintf(fid, ', %17.15g ', t(n));
end
fprintf(fid, '};\n');

%% interp to quadrature
Nq = length(rq);

fprintf(fid,'int p_Nq = %d;\n',Nq);

fprintf(fid, 'dfloat p_wq[%d] = {', Nq);
fprintf(fid, '%17.15g ', w(1));
for n=2:Nq
    fprintf(fid, ', %17.15g ', w(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_rq[%d] = {', Nq);
fprintf(fid, '%17.15g ', rq(1));
for n=2:Nq
    fprintf(fid, ', %17.15g ', rq(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_sq[%d] = {', Nq);
fprintf(fid, '%17.15g ', sq(1));
for n=2:Nq
    fprintf(fid, ', %17.15g ', sq(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_tq[%d] = {', Nq);
fprintf(fid, '%17.15g ', tq(1));
for n=2:Nq
    fprintf(fid, ', %17.15g ', tq(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_Vq[%d][%d] = {', Nq, Np);
for m=1:Nq
    fprintf(fid, '{%17.15g ', Vq(m,1));
    for n=2:Np
        fprintf(fid, ', %17.15g ', Vq(m,n));
    end
    if(m<Nq)
        fprintf(fid, '},\n');
    else
        fprintf(fid, '}};\n');
    end
end

%% triangle quadrature

[rfq sfq wfq] = Cubature2D(2*N);

Nfq = length(rfq);
fprintf(fid,'int p_Nfq = %d;\n',Nfq);

fprintf(fid, 'dfloat p_wfq[%d] = {', Nfq);
fprintf(fid, '%17.15g ', wfq(1));
for n=2:Nfq
    fprintf(fid, ', %17.15g ', wfq(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_rfq[%d] = {', Nfq);
fprintf(fid, '%17.15g ', rfq(1));
for n=2:Nfq
    fprintf(fid, ', %17.15g ', rfq(n));
end
fprintf(fid, '};\n');

fprintf(fid, 'dfloat p_sfq[%d] = {', Nfq);
fprintf(fid, '%17.15g ', sfq(1));
for n=2:Nfq
    fprintf(fid, ', %17.15g ', sfq(n));
end
fprintf(fid, '};\n');

%%

% fprintf(fid, 'dfloat p_Dr[%d][%d] = {', Np, Np);
% for m=1:Np
%     fprintf(fid, '{%17.15g ', Dr(m,1));
%     for n=2:Np
%         fprintf(fid, ', %17.15g ', Dr(m,n));
%     end
%     if(m<Np)
%         fprintf(fid, '},\n');
%     else
%         fprintf(fid, '}};\n');
%     end
% end
% 
% fprintf(fid, 'dfloat p_Ds[%d][%d] = {', Np, Np);
% for m=1:Np
%     fprintf(fid, '{%17.15g ', Ds(m,1));
%     for n=2:Np
%         fprintf(fid, ', %17.15g ', Ds(m,n));
%     end
%     if(m<Np)
%         fprintf(fid, '},\n');
%     else
%         fprintf(fid, '}};\n');
%     end
% end
% 
% fprintf(fid, 'dfloat p_Dt[%d][%d] = {', Np, Np);
% for m=1:Np
%     fprintf(fid, '{%17.15g ', Dt(m,1));
%     for n=2:Np
%         fprintf(fid, ', %17.15g ', Dt(m,n));
%     end
%     if(m<Np)
%         fprintf(fid, '},\n');
%     else
%         fprintf(fid, '}};\n');
%     end
% end
% 
% %fprintf(fid, '// blaaaaah\n');
% fprintf(fid, 'dfloat p_LIFT[%d][%d] = {', Np, Nfp*Nfaces);
% for m=1:Np
%     fprintf(fid, '{%17.15g ', LIFT(m,1));
%     for n=2:Nfp*Nfaces
%         fprintf(fid, ', %17.15g ', LIFT(m,n));
%     end
%     if(m<Np)
%         fprintf(fid, '},\n');
%     else
%         fprintf(fid, '}};\n');
%     end
% end

fprintf(fid, 'int p_Fmask[%d][%d] = {', Nfaces, Nfp);
for m=1:Nfaces
    fprintf(fid, '{%d ', Fmask(1,m));
    for n=2:Nfp
        fprintf(fid, ', %d ', Fmask(n,m));
    end
    if(m<Nfaces)
        fprintf(fid, '},\n');
    else
        fprintf(fid, '}};\n');
    end
end

%% bern stuff
[VB Vr Vs Vt V1 V2 V3 V4 id] = bern_basis_tet(N,r,s,t);
  
% V, invV (convert from nodal to Bern basis)
write_matrix(fid,'VB',VB)
write_matrix(fid,'invVB',inv(VB))

DrB = VB\Vr; DrB(abs(DrB)<1e-8) = 0;
DsB = VB\Vs; DsB(abs(DsB)<1e-8) = 0;
DtB = VB\Vt; DtB(abs(DtB)<1e-8) = 0;

% DrstB matices for convenience
write_matrix(fid,'DrB',DrB)
write_matrix(fid,'DsB',DsB)
write_matrix(fid,'DtB',DtB)

D1 = VB\V1;  D1(abs(D1)<1e-8) = 0;
D2 = VB\V2;  D2(abs(D2)<1e-8) = 0;
D3 = VB\V3;  D3(abs(D3)<1e-8) = 0;
D4 = VB\V4;  D4(abs(D4)<1e-8) = 0;
VBq = bern_basis_tet(N,rq,sq,tq); 
M = VBq'*diag(w)*VBq; %invM = inv(M);

% D1 to D4
D1_ids = zeros(size(D1,1),4); D1_vals = zeros(size(D1,1),4);
D2_ids = zeros(size(D2,1),4); D2_vals = zeros(size(D2,1),4);
D3_ids = zeros(size(D3,1),4); D3_vals = zeros(size(D3,1),4);
D4_ids = zeros(size(D4,1),4); D4_vals = zeros(size(D4,1),4);
for i = 1:size(D1,1)
    rowids = find(D1(i,:)); extraids = setdiff(1:size(D1,2),rowids);
    D1_ids(i,:) = [rowids extraids(1:4-length(rowids))];    
    D1_vals(i,:) = D1(i,D1_ids(i,:));
    
    rowids = find(D2(i,:)); extraids = setdiff(1:size(D2,2),rowids);
    D2_ids(i,:) = [rowids extraids(1:4-length(rowids))];    
    D2_vals(i,:) = D2(i,D2_ids(i,:));

    rowids = find(D3(i,:)); extraids = setdiff(1:size(D3,2),rowids);
    D3_ids(i,:) = [rowids extraids(1:4-length(rowids))];    
    D3_vals(i,:) = D3(i,D3_ids(i,:));

    rowids = find(D4(i,:)); extraids = setdiff(1:size(D4,2),rowids);
    D4_ids(i,:) = [rowids extraids(1:4-length(rowids))];    
    D4_vals(i,:) = D4(i,D4_ids(i,:));
end

diff = norm(D1_vals-D2_vals,'fro') + norm(D1_vals-D3_vals,'fro') + norm(D1_vals-D4_vals,'fro');
if diff>1e-8
    keyboard
end

% all values are the same! 
write_int_matrix(fid,'D1_ids',D1_ids-1); write_matrix(fid,'D_vals',D1_vals)
write_int_matrix(fid,'D2_ids',D2_ids-1); %write_matrix(fid,'D2_vals',D2_vals)
write_int_matrix(fid,'D3_ids',D3_ids-1); %write_matrix(fid,'D3_vals',D3_vals)
write_int_matrix(fid,'D4_ids',D4_ids-1); %write_matrix(fid,'D4_vals',D4_vals)

% compressed ID storage
sk = 1;
for i = 0:N
    for j = 0:N-i
        for k = 0:N-i-j
            l = N-i-j-k;
            idi(sk) = i; idj(sk) = j; idk(sk) = k; idl(sk) = l;
            sk = sk + 1;
        end
    end
end
idi = idi(id); idj = idj(id); idk = idk(id); idl = idl(id);
ids = [idi;idj;idk;idl];

write_int_matrix(fid,'D_iids',idi);
write_int_matrix(fid,'D_jids',idj);
write_int_matrix(fid,'D_kids',idk);

% LIFT, LIFTrst 
Emat = zeros(Np, Nfaces*Nfp);
Mf = zeros(Np, Nfp);
for face=1:Nfaces
    % WARNING: ASSUMES ORDERING OF FACE NODES/FMASK    
    [rfq sfq wf] = Cubature2D(2*N);         
    Vfq = bern_basis_tri(N,rfq,sfq);
    massFace = Vfq'*diag(wf)*Vfq;
    
    idr = Fmask(:,face)+1; % undoes previous 0 indexing...
    idc = (face-1)*Nfp+1:face*Nfp;
    
    Emat(idr, idc) = Emat(idr, idc) + massFace;    
    if face==1
        Mf(idr,idc)=massFace;
    end
end
LIFTB = M\Emat; LIFTB(abs(LIFTB)<1e-8) = 0;
write_matrix(fid,'LIFTB',LIFTB)

%% compressed versions of LIFT matrices: store L0 + EEL reduction matrix (sparse)

% decompose lift operators
LIFTf = M\Mf; LIFTf(abs(LIFTf)<1e-8) = 0;
LIFTB(abs(LIFTB)<1e-8) = 0;

% get permutations of rows for 2nd-4th cols of LIFT
p = zeros(Np,Nfaces-1);
for f = 1:Nfaces-1
    ids = (1:Nfp) + f*Nfp;
    for i = 1:Np % match rows of first cols with second
        diff = kron(ones(Np,1),LIFTf(i,:)) - LIFTB(:,ids);
        [~,j] = min(sum(abs(diff),2));
        p(j,f) = i;
    end
end

[r2D s2D] = Nodes2D(N); [r2D s2D] = xytors(r2D,s2D);
Vtri = bern_basis_tri(N,r2D,s2D);

% 2d degree ops
ED = {}; 
for i = 0:N     
    EDi = bern_basis_tri(N,r2D,s2D)\bern_basis_tri(N-i,r2D,s2D);
    EDi(abs(EDi)<1e-8) = 0;
    ED{i+1} = EDi;    
end
EE = blkdiag(ED{:}); EE = EE';

% build in scalings
EDL = cell(N+1,1);
for i = 0:N
    l(i+1) = (-1)^i*nchoosek(N,i)/(1+i); % coefficients in the lift decomposition
    EDL{i+1} = ED{i+1}*l(i+1);
end

% lift reduction operator
EL = [];
for i = 0:N
    EL = [EL;EDL{i+1}'];
end
EEL=([EL EL(p(:,1),:) EL(p(:,2),:) EL(p(:,3),:)]);

% scaled degree elevation matrix for all faces
EEL(abs(EEL)<1e-6) = 0;
EEL_nnz = max(sum(abs(EEL)>0,2)); % max nnz in a row for EEL - should be Nfp + 2

EEL_ids = zeros(size(EEL,1),EEL_nnz); 
EEL_vals = zeros(size(EEL,1),EEL_nnz);
for i = 1:size(EEL,1)
    rowids = find(EEL(i,:)); extraids = setdiff(1:size(EEL,2),rowids);
    EEL_ids(i,:) = [rowids extraids(1:EEL_nnz-length(rowids))];    
    EEL_vals(i,:) = EEL(i,EEL_ids(i,:));
end

fprintf(fid, 'int p_EEL_nnz = %d;\n', EEL_nnz);
write_matrix(fid,'EEL_vals',EEL_vals)
write_int_matrix(fid,'EEL_ids',EEL_ids-1); % decrement to 0 index

% % compute L0
% [rplus splus] = Nodes2D(N+1); [rplus splus] = xytors(rplus,splus);
% Eplus = bern_basis_tri(N+1,rplus,splus)\bern_basis_tri(N,rplus,splus);
% L0 = Eplus'*Eplus * (N+1)^2/2; 
L0 = LIFTB(1:Nfp,1:Nfp);

% extract sparsity pattern
L0_nnz = min(Nfp,7); 
L0_ids = zeros(size(L0,1),L0_nnz);
L0_vals = zeros(size(L0,1),L0_nnz);
for i = 1:size(L0,1)
    rowids = find(L0(i,:)); extraids = setdiff(1:size(L0,2),rowids);
    L0_ids(i,:) = [rowids extraids(1:L0_nnz-length(rowids))];    
    L0_vals(i,:) = L0(i,L0_ids(i,:));
end

write_matrix(fid,'L0_vals',L0_vals)
write_int_matrix(fid,'L0_ids',L0_ids-1); % decrement to 0 index

for i = 0:N
    cEL(i+1) = (-1)^i*nchoosek(N,i)/(1+i);
end
% fprintf(fid, 'dfloat p_cEL[%d] = { ', N+1);
% fprintf(fid, '%17.15g ', cEL(1));
% for n=2:N+1
%     fprintf(fid, ', %17.15g ', cEL(n));
% end
% fprintf(fid, '};\n');


%% close file

fprintf(fid, sprintf('\n#endif\n'));
fclose(fid)


function write_int_matrix(fid,name,A)

[rows cols] = size(A);
fprintf(fid, 'int p_%s[%d][%d] = {', name,rows,cols);
for m=1:rows
    fprintf(fid, '{%d ', A(m,1));
    for n=2:cols
        fprintf(fid, ', %d ', A(m,n));
    end
    if(m<rows)
        fprintf(fid, '},\n');
    else
        fprintf(fid, '}};\n');
    end
end


function write_matrix(fid,name,A)

[rows cols] = size(A);
fprintf(fid, 'dfloat p_%s[%d][%d] = {', name,rows,cols);
for m=1:rows
    fprintf(fid, '{%17.15g ', A(m,1));
    for n=2:cols
        fprintf(fid, ', %17.15g ', A(m,n));
    end
    if(m<rows)
        fprintf(fid, '},\n');
    else
        fprintf(fid, '}};\n');
    end
end
