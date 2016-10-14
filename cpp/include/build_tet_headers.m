% writes header files for nodal bases

function build_tet_headers

% for N = 1:7
%     write_headers(N);
% end

% single file containing all node data
fid = fopen(sprintf('NodeData.h'), 'w');
fprintf(fid, sprintf('#ifndef NODEDATA\n\n'));
fprintf(fid, sprintf('#define NODEDATA 1\n\n'));
for N = 1:21 % to 3*N for N = 7
    fprintf('on N = %d\n',N)
    
%     [r1D00 wq1D00] = JacobiGQ(0,0,N);
%     [r1D10 wq1D10] = JacobiGQ(1,0,N);
%     [r1D20 wq1D20] = JacobiGQ(2,0,N);
%     write_vector(fid,sprintf('r1D00_N%d',N),r1D00);
%     write_vector(fid,sprintf('r1D10_N%d',N),r1D10);
%     write_vector(fid,sprintf('r1D20_N%d',N),r1D20);    
%     write_vector(fid,sprintf('wq1D00_N%d',N),wq1D00);
%     write_vector(fid,sprintf('wq1D10_N%d',N),wq1D10);
%     write_vector(fid,sprintf('wq1D20_N%d',N),wq1D20);        
    
    [r s t] = Nodes3D(N); [r s t] = xyztorst(r,s,t);
    write_vector(fid,sprintf('r_N%d',N),r);
    write_vector(fid,sprintf('s_N%d',N),s);
    write_vector(fid,sprintf('t_N%d',N),t);    
    
    [rq sq tq wq] = tet_cubature(N);
    fprintf(fid,'int p_Nq_N%d = %d;\n',N,length(wq));
    write_vector(fid,sprintf('rq_N%d',N),rq);
    write_vector(fid,sprintf('sq_N%d',N),sq);
    write_vector(fid,sprintf('tq_N%d',N),tq);
    write_vector(fid,sprintf('wq_N%d',N),wq);
    
    [rfq sfq wfq] = Cubature2D(N);
    fprintf(fid,'int p_Nfq_N%d = %d;\n',N,length(wfq));
    write_vector(fid,sprintf('rfq_N%d',N),rfq);
    write_vector(fid,sprintf('sfq_N%d',N),sfq);
    write_vector(fid,sprintf('wfq_N%d',N),wfq);
end
% close file
fprintf(fid, sprintf('\n#endif\n'));
fclose(fid)


function write_headers(in_N)

if nargin==0
    in_N=1;
end

Globals3D;

NODETOL = 1e-8;

N = in_N;

Np = (N+1)*(N+2)*(N+3)/6;
Nfp = (N+1)*(N+2)/2;
Nfaces = 4;

[x,y,z] = Nodes3D(N);
[r,s,t] = xyztorst(x,y,z);

% find all the nodes that lie on each edge
fmask1   = find( abs(t+1) < NODETOL)';
fmask2   = find( abs(s+1) < NODETOL)';
fmask3   = find( abs(r+s+t+1) < NODETOL)';
fmask4   = find( abs(r+1) < NODETOL)';
Fmask  = [fmask1;fmask2;fmask3;fmask4]';

% interp to quadrature points
[rq sq tq wq] = tet_cubature(2*N+1);

% correct to 0-index
Fmask = Fmask-1;

%%

fid = fopen(sprintf('data3dN0%d.h', N), 'w');

fprintf(fid, sprintf('#ifndef DATA3dN0%d \n\n', N));
fprintf(fid, sprintf('#define DATA3dN0%d 1\n\n', N));

write_vector(fid,'r',r);
write_vector(fid,'s',s);
write_vector(fid,'t',t);

% quadrature
Nq = length(rq);
fprintf(fid,'int p_Nq = %d;\n',Nq);
write_vector(fid,'wq',wq);
write_vector(fid,'rq',rq);
write_vector(fid,'sq',sq);
write_vector(fid,'tq',tq);

% face quadrature
[rfq sfq wfq] = Cubature2D(2*N);
Nfq = length(rfq);
fprintf(fid,'int p_Nfq = %d;\n',Nfq);
write_vector(fid,'wfq',wfq);
write_vector(fid,'rfq',rfq);
write_vector(fid,'sfq',sfq);
write_int_matrix(fid,'Fmask',Fmask');

% close file
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
fprintf(fid, 'double p_%s[%d][%d] = {', name,rows,cols);
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

function write_vector(fid,name,A)

rows = length(A(:));
fprintf(fid, 'double p_%s[%d] = {', name,rows);
fprintf(fid, '%17.15g ', A(1));
for m=2:rows
    fprintf(fid, ', %17.15g ', A(m));
end
fprintf(fid, '};\n');

