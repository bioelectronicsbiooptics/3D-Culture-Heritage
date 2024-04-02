% STL_Decode_v8_SYJ
%% Fastq text scan Jikji6
% SYJ Edited
% filename = input("filename?: ",'s');
filename ='Jikji6' % 'Jikji6' 'seokga' 'pisa' 'liberty' 'tripitaka' 'venus' 'goldengate'

tStart = tic;
tic 
for i = 1 : 2
    % fid = fopen(filename+"_"+i+".fastq");
    fid = fopen("6_S2_L001_R"+i+"_001"+".fastq");
    % fid = fopen("G_"+i+".fastq");

    warning('off','all')
    data = textscan(fid,'%s','Delimiter','\n','CommentStyle',{'+','@'});
    warning('on','all')
    if i == 1
        seq_1 = data{1}(2:2:end);
    else
        seq_2 = data{1}(2:2:end);
    end
    fclose(fid);
end
seq = [seq_1;seq_2];

addpath('Fastq Function2')
T(1)= toc; 
"Paired End Sequence Read: " + T(1)
%% Evert Sequence
tic
ForwardSection = cell(length(seq),1);

chunkSize = 1000;
numChunks = ceil(length(seq) / chunkSize);
ForwardSection = cell(size(seq));

for i = 1:numChunks
    startIndex = (i - 1) * chunkSize + 1;
    endIndex = min(i * chunkSize, length(seq));

    for j = startIndex:endIndex
        ForwardSection{j} = seq{j}(1:21);
    end
end

[F_primers,R_primers] = PrimerSelect(filename);
[tevertSeq, ReversedSeqLog] = ReverseEvert(seq,ForwardSection,R_primers,F_primers);
seq(logical(ReversedSeqLog)) = tevertSeq;
T(2)= toc; 
"Reverse Complemantary Sequence: "+ T(2)
%% Primer Sorting
tic
PrimerSortingSeq = PrimerSort(seq,F_primers,R_primers);
clear tevertSeq ReversedSeqLog ForwardSection
T(3)=toc;
"Primer Sorting: "+T(3)
%% Key finding
tic
[KeyMap,KeymatDNA] = KeyFind(PrimerSortingSeq);
k = 156;
n = k+8;
m = length(de2bi(n));
[OligoLen,idxlen,ncolpo,nbitpo,ncolco,nbitco] = Keysort(KeyMap,KeymatDNA,k,filename);
T(4)= toc; 
"Key finding: " + T(4)
%% Consensus Sequencing
tic
DNAmatrix = cell2mat(PrimerSortingSeq);
MIndex = Consensus2(DNAmatrix,idxlen,OligoLen);

SeqtoBin = Binarization(MIndex);
T(5)=toc;
"Consensus Sequencing & binerization: " +T(5)

Depth = zeros(length(MIndex),1);
for i = 1 : length(MIndex)
    Depth(i) = size(MIndex{i},1);
end
"Depth mean: "+mean(Depth)
%% RS decoding
tic
[rsECC,RSerrinf] = RSdecoding(SeqtoBin,k,m,n);
GRS_err_indices = find(RSerrinf<0);
rsECC(RSerrinf<0,:) = zeros(length(find(RSerrinf<0)),size(rsECC,2));
rsECC(RSerrinf>4,:) = zeros(length(find(RSerrinf>4)),size(rsECC,2));
for a = 1 : OligoLen
    if ~isempty(find(rsECC(a,:)>1, 1))
        rsECC(a,:) = zeros(1,size(rsECC,2));
    end
end
payload = rsECC(:,1:k-idxlen);

if strcmp(filename,'pisa') || strcmp(filename,'seokga') || strcmp(filename,'liberty')
    Temp_mat = payload(1:end-2,:);
else
Temp_mat = payload(1:end-1,:);
end

T(6)=toc;
"RS decoding: " +T(6)
%% Null finding
tic
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

[data_Ndecode,xor_Ndecode,dNc,xNc] = NullFind(seq,Null_count,Null_xor_count,s,F_primers);
d = {'00','01','10','11'};
[~,x] = ismember(data_Ndecode(21:end),'ATGC');
chr_out = cell2mat(d(x));
bindataN = chr_out-'0';

[~,y] = ismember(xor_Ndecode(21:end),'ATGC');
chr_out = cell2mat(d(y));
binxorN = chr_out-'0';

Depth(floor((numel(Temp_mat)*2/3)/s)) = size(cell2mat(dNc'),1);
Depth(end-1) = size(cell2mat(xNc'),1);
T(7) = toc;
"Null finding " + T(7)


fid = fopen(filename+".txt");
[scan,count] = fscanf(fid,'%s');
orgseq = reshape(scan,[length(scan)/count count])';

if sum(data_Ndecode(21:end) == orgseq(floor((numel(Temp_mat)*2/3)/s),21:length(data_Ndecode))) == length(data_Ndecode)-20
    "Data part Null Decoding Complete"
end

if sum(xor_Ndecode(21:end) == orgseq(end-1,21:length(xor_Ndecode))) == length(xor_Ndecode)-20
    "XOR part Null Decoding Complete"
end

%% XOR
tic

matsz = k-idxlen;
Temp_mat(floor((numel(Temp_mat)*2/3)/s),1:length(bindataN)) = bindataN; 
Temp_mat(end,1:length(binxorN)) = binxorN;
if size(Depth,1) > size(Temp_mat,1)
    Depth(end) = [];
end
T(8)= toc;
"xxxxxx " + T(8)
%%
tic
try
    ECC_tri_part = Temp_mat(1:floor(numel(Temp_mat)*2/3/s),:);
    ECC_payload = reshape(ECC_tri_part',[numel(ECC_tri_part) 1])';
    ECC_payload(end-Null_count+1:end) = [];
    [tc,tp] = tridecode(ECC_payload,ncolpo,nbitpo,ncolco,nbitco,filename);
    trid = trierrXOR(Temp_mat,tc,tp,Null_count,Null_xor_count,ncolpo,nbitpo,ncolco,nbitco);
    trid = triangulation(tc,tp);
    figure
    trimesh(trid)
    title('Global RS')
catch

end

Depthcut = 16;
ErrD =find(Depth<=Depthcut);

try
    Temp_mat(ErrD(RSerrinf(ErrD)~=0),:) = zeros(numel(ErrD(RSerrinf(ErrD)~=0)),matsz);
    [ECC_tri_part,~,~,~,~,~,~] = ...
        XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);

    ECC_payload = reshape(ECC_tri_part',[numel(ECC_tri_part) 1])';
    ECC_payload(end-Null_count+1:end) = [];
    [tc,tp] = tridecode(ECC_payload,ncolpo,nbitpo,ncolco,nbitco,filename);
    trid = trierrXOR(Temp_mat,tc,tp,Null_count,Null_xor_count,ncolpo,nbitpo,ncolco,nbitco);
    figure
    trimesh(trid)
    title('Depth cut')
catch

end

try
    Temp_mat(ErrD(RSerrinf(ErrD)~=0),:) = zeros(numel(ErrD(RSerrinf(ErrD)~=0)),matsz);
    [ECC_tri_part,~,tri_part,xor_part,XOR_ECC_A,XOR_ECC_B,XOR_ECC_C] = ...
    XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);

    [XORMIndex,Del_xor] = nestedXOR(MIndex,ECC_tri_part,ErrD,RSerrinf,XOR_ECC_A, ...
        XOR_ECC_B,XOR_ECC_C,payload,tri_part,xor_part,matsz,Null_count);

    Temp_mat = reRSdec(XORMIndex,Temp_mat,k,m,n,Del_xor,matsz);
    [~,missing_data,~,~,~,~,~] = ...
        XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);

    XORmData = MIndex(missing_data);
    Temp_mat = reRSdec(XORmData,Temp_mat,k,m,n,missing_data,matsz);

    [ECC_tri_part,Local_RS_final_Index,~,~,~,~,~] = ...
        XORrun(Temp_mat,s,Null_count,Null_xor_count,ncolpo,ncolco,nbitpo,nbitco);
    ECC_payload = reshape(ECC_tri_part',[numel(ECC_tri_part) 1])';
    ECC_payload(end-Null_count+1:end) = [];
    [tc,tp] = tridecode(ECC_payload,ncolpo,nbitpo,ncolco,nbitco,filename);
    trid = trierrXOR(Temp_mat,tc,tp,Null_count,Null_xor_count,ncolpo,nbitpo,ncolco,nbitco);
    length(Del_xor),length(missing_data),length(Local_RS_final_Index)
    figure
    trimesh(trid)
    title('Local RS')
catch
    'Decoding Failed'
end
T(9)= toc;
"tri_tri " + T(9)
tMul = sum(T)
tMul/60