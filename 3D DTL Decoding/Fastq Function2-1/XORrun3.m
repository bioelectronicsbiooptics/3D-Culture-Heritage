function [xtp,md,tp,xp,XA,XB,XC] = XORrun3(matrix,s,Nc,Nxc,ncp,ncc,nbp,nbc)
%trid = trierrXOR(Temp_mat,tc,tp,Null_count,Null_xor_count,ncolpo,nbitpo,ncolco,nbitco); 
% matrix =Temp_mat;
% Nc = Null_count;
% Nxc =Null_xor_count;
% ncp = ncolpo;
% ncc = nbitpo;
% nbp= ncolco;
% nbc = nbitco;

tp = matrix(1:floor(numel(matrix)*2/3/s),:);
xp = matrix(floor(numel(matrix)*2/3/s)+1:end-1,:);%-1 check!!!!

row_part = reshape(tp',[numel(tp) 1])';
row_part(end-Nc+1:end) = [];

if rem(ncp*3*nbp+ncc*3*nbc,2) ~= 0
    row_part(end) = [];
end

XC = reshape(xp',[numel(xp) 1])';
XC(end-Nxc+1:end) = [];

XA = row_part(1:ceil(length(row_part)/2));
XB = row_part(length(XA)+1:end);
if rem(length(row_part),2) ~= 0
    XB(end+1) = 0;
end
XOR_XA_XB = [xor(XB,XC) xor(XA,XC)];
XOR_XC = xor(XA,XB);

ECC_example = [XOR_XA_XB zeros(1,Nc)];
ECC_ex_matrix = reshape(ECC_example,[s numel(ECC_example)/s])';

ECC_example2 = [XOR_XC zeros(1,Nxc)];
ECC_ex_matrix2 = reshape(ECC_example2,[s numel(ECC_example2)/s])';
ECC_ex_matrix3 = [ECC_ex_matrix ; ECC_ex_matrix2];
md_t = find(sum(tp')==0); % RS가 fail 인걸해야 (위에서 선언했음)
md_x = find(sum(xp')==0); % RS가 fail 인걸해야 (위에서 선언했음)
md = find(sum(matrix')==0); % RS가 fail 인걸해야 (위에서 선언했음)
%치환되는 부분!
tp(md_t,:) = ECC_ex_matrix(md_t,:);
xp(md_x,:) = ECC_ex_matrix2(md_x,:);
matrix(md,:) = ECC_ex_matrix3(md,:);
xtp = matrix;
