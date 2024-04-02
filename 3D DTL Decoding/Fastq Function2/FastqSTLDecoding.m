function trid = FastqSTLDecoding(filename)

seq = fastqTextscan(filename);

ForwardSection = cell(length(seq),1);
parfor i = 1 : length(seq)
    ForwardSection{i} = seq{i}(1:21);
end
[F_primers,R_primers] = PrimerSelect(filename);
[tevertSeq, ReversedSeqLog] = ReverseEvert(seq,ForwardSection,R_primers);
seq(logical(ReversedSeqLog)) = tevertSeq;

PrimerSortingSeq = PrimerSort(seq,F_primers,R_primers);

[KeyMap,KeymatDNA] = KeyFind(PrimerSortingSeq);
k = 156;
n = k+8;
m = length(de2bi(n));
[OligoLen,idxlen,ncolpo,nbitpo,ncolco,nbitco] = Keysort(KeyMap,KeymatDNA,k,filename);

DNAmatrix = cell2mat(PrimerSortingSeq);
MIndex = Consensus(DNAmatrix,idxlen,OligoLen);

Depth = zeros(length(MIndex),1);
for i = 1 : length(MIndex)
    Depth(i) = size(MIndex{i},1);
end

if strcmp(filename,'pisa') || strcmp(filename,'seokga') || strcmp(filename,'liberty')
    [MIndex, Map] = MapLocation(MIndex);
end    

SeqtoBin = Binarization(MIndex);
[rsECC,RSerrinf] = RSdecoding(SeqtoBin,k,m,n);
rsECC(RSerrinf<0,:) = zeros(length(find(RSerrinf<0)),size(rsECC,2));
rsECC(RSerrinf>4,:) = zeros(length(find(RSerrinf>4)),size(rsECC,2));
for a = 1 : OligoLen
    if ~isempty(find(rsECC(a,:)>1, 1))
        rsECC(a,:) = zeros(1,size(rsECC,2));
    end
end
payload = rsECC(:,1:k-idxlen);
Temp_mat = payload(1:end-1,:);

if rem(ncolpo*3*nbitpo+ncolco*3*nbitco,2) == 1
    Nullkey_1 = ncolpo*3*nbitpo+ncolco*3*nbitco+1;
else
    Nullkey_1 = ncolpo*3*nbitpo+ncolco*3*nbitco;
end
s = size(Temp_mat,2);

Nullkey_2 = Nullkey_1/2;
Null_count = 0;
while rem(Nullkey_1,s) ~= 0
    Nullkey_1 = Nullkey_1 + 1;
    Null_count = Null_count + 1;
end
Null_xor_count = 0;
while rem(Nullkey_2,s) ~= 0
    Nullkey_2 = Nullkey_2 + 1;
    Null_xor_count = Null_xor_count + 1;
end
[data_Ndecode,xor_Ndecode] = NullFind(seq,Null_count,Null_xor_count,s,F_primers);

d = {'00','01','10','11'};
[~,x] = ismember(data_Ndecode(21:end),'ATGC');
chr_out = cell2mat(d(x));
bindataN = chr_out-'0';

[~,y] = ismember(xor_Ndecode(21:end),'ATGC');
chr_out = cell2mat(d(y));
binxorN = chr_out-'0';

Temp_mat(floor((numel(Temp_mat)*2/3)/s),1:length(bindataN)) = bindataN; 
Temp_mat(end,1:length(binxorN)) = binxorN;
Temp_mat(Depth<=16,:) = zeros(numel(find(Depth<=16)),140);

ErrD =find(Depth<=16);
ErrD(RSerrinf(ErrD)~=0);
Temp_mat(ErrD(RSerrinf(ErrD)~=0),:) = zeros(numel(ErrD(RSerrinf(ErrD)~=0)),140);

[ECC_tri_part,~,tri_part,xor_part,XOR_ECC_A,XOR_ECC_B,XOR_ECC_C] = ...
    XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);
matsz = k-idxlen;
if ~isempty(ErrD)
    [XORMIndex,Del_xor] = nestedXOR(MIndex,ECC_tri_part,ErrD,RSerrinf,XOR_ECC_A, ...
        XOR_ECC_B,XOR_ECC_C,payload,tri_part,xor_part,matsz,Null_count);
    
    Temp_mat = reRSdec(XORMIndex,Temp_mat,k,m,n,Del_xor);
    [~,missing_data,~,~,~,~,~] = ...
        XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);
    
    XORmData = MIndex(missing_data);
    Temp_mat = reRSdec(XORmData,Temp_mat,k,m,n,missing_data);
end
[ECC_tri_part,~,~,~,~,~,~] = ...
    XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);

ECC_payload = reshape(ECC_tri_part',[numel(ECC_tri_part) 1])';
ECC_payload(end-Null_count+1:end) = [];

trid = tridecode(ECC_payload,ncolpo,nbitpo,ncolco,nbitco,filename);
trimesh(trid)