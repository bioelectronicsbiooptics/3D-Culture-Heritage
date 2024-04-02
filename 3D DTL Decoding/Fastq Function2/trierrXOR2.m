function trid = trierrXOR(mat,tc,tp,Ncount,Nxcount,ncp,nbp,ncc,nbc)
s = size(mat,2);
tri_part = mat(1:floor(numel(mat)*2/3/s),:);
xor_part = mat(floor(numel(mat)*2/3/s)+1:end,:);

row_tri_part = reshape(tri_part',[numel(tri_part) 1])';
row_tri_part(end-Ncount+1:end) = [];

if rem(ncp*3*nbp+ncc*3*nbc,2) ~= 0
    row_tri_part(end) = [];
end

XOR_ECC_C = reshape(xor_part',[numel(xor_part) 1])';
XOR_ECC_C(end-Nxcount+1:end) = [];

XOR_ECC_A = row_tri_part(1:ceil(length(row_tri_part)/2));
XOR_ECC_B = row_tri_part(length(XOR_ECC_A)+1:end);
if rem(length(row_tri_part),2) ~= 0
    XOR_ECC_B(end+1) = 0;
end
XOR_ECC_matrix = [xor(XOR_ECC_B,XOR_ECC_C) xor(XOR_ECC_A,XOR_ECC_C)];
%%
XOR_ECC_C_Matrix = xor(XOR_ECC_A,XOR_ECC_B);
%%         
%% resize

if numel(find(~tp')) > 0
    XOR_rep = find(~tp');
    XOR_prng = [XOR_rep*nbp-nbp+1 XOR_rep*nbp];
    eccxorpo = zeros(size(XOR_prng,1),nbp);
    for a = 1 : size(XOR_prng,1)
        eccxorpo(a,:) = XOR_ECC_matrix(XOR_prng(a,1) : XOR_prng(a,2));
    end
    respo = bi2de(eccxorpo)/10^16;
    transposePO = tp';
    transposePO(XOR_rep) = respo;
    tp = transposePO';
end

if numel(find(~tc')) > 0 || (numel(find((tc'>size(tp,1)))))
    XOR_req = [find(~tc'); find(tc'>size(tp,1))];
    XOR_rng = [XOR_req*nbc-nbc+1 XOR_req*nbc];
    coXOR1 = XOR_ECC_C;
    coXOR2 = XOR_ECC_A;
    if rem(length(row_tri_part),2) ~= 0
        coXOR1(end) = [];
        coXOR2(end) = [];
    end
    XORrest = xor(coXOR1,coXOR2);
    XR = [];
    for xr = 1 : size(XOR_rng,1)
        XC = XORrest(XOR_rng(xr,1):XOR_rng(xr,2));
        XR = [XR; XC];
    end
    transposeCO = tc';
    transposeCO(XOR_req) = bi2de(XR);
    tc = transposeCO';

    nullsz = [find(tc > ncp); find(tc == 0)];
    for a = 1 : length(nullsz)
        while nullsz(a) > ncc
            nullsz(a) = nullsz(a) -ncc;
        end
    end
    size(tc(unique(nullsz),:),1)
    tc(unique(nullsz),:) = [];
end

trid = triangulation(tc,tp);