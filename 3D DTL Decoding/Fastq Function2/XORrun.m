function [xtp,md,tp,xp,XA,XB,XC] = XORrun(matrix,s,Nc,Nxc,ncp,ncc,nbp,nbc)
tp = matrix(1:floor(numel(matrix)*2/3/s),:);
xp = matrix(floor(numel(matrix)*2/3/s)+1:end,:);

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
XOR_matrix = [xor(XB,XC) xor(XA,XC)];
XOR_XC_matrix = xor(XA,XB);

ECC_example = [XOR_matrix zeros(1,Nc)];
ECC_ex_matrix = reshape(ECC_example,[s numel(ECC_example)/s])';
md = find(sum(tp')==0);
tp(md,:) = ECC_ex_matrix(md,:);
xtp = tp;