function CellSeq = Consensus(DNAmatrix, idxlen, OligoLen)
d = {'00','01','10','11'};
IndexDNA = DNAmatrix(:,98-ceil(idxlen/2)+1:98);
Indexbin = zeros(length(IndexDNA),size(IndexDNA,2)*2);
for r = 1 : size(IndexDNA,1)
    if (strfind(IndexDNA(r,:),'N'))>0
        IndexDNA(r,:) = repmat('A',1,size(IndexDNA,2));
        continue
    end
    [~,x] = ismember(IndexDNA(r,:),'ATGC');
    chr_out = cell2mat(d(x));
    Indexbin(r,:) = chr_out-'0';
end
if rem(length(de2bi(OligoLen)),2) ~= 0
    if rem(size(Indexbin,2),2) == 0
        Indexbin = Indexbin(:,2:end);
    end
end
Indexdec = bi2de(Indexbin);

CellSeq = cell(OligoLen,1);
if OligoLen>=1000
    for i = 1 : floor(OligoLen/1000)
        for a = 1000*(i-1)+1 : 1000*i
            CellSeq{a} = DNAmatrix((Indexdec == a),:);
        end
    end
    for j = a(end)+1 : OligoLen
        CellSeq{j} = DNAmatrix((Indexdec == j),:);
    end
else
    for a = 1 : OligoLen
        CellSeq{a} = DNAmatrix((Indexdec == a),:);
    end
end