function CellSeq = Consensus(DNAmatrix, idxlen, OligoLen)

d = {'00','01','10','11'};
IndexDNA = DNAmatrix(:,98-ceil(idxlen/2)+1:98);
Indexbin = zeros(length(IndexDNA),size(IndexDNA,2)*2);
for r = 1:size(IndexDNA, 1)
    if any(IndexDNA(r,:) == 'N')
         IndexDNA(r,:) = repmat('A', 1, size(IndexDNA, 2));
    end
end
parfor r = 1:size(IndexDNA, 1)
         [~, x] = ismember(IndexDNA(r,:), 'ATGC');
        chr_out = cell2mat(d(x));
        Indexbin(r,:) = chr_out - '0';
end
if rem(length(de2bi(OligoLen)),2) ~= 0
    if rem(size(Indexbin,2),2) == 0
        Indexbin = Indexbin(:,2:end);
    end
end
Indexdec = bi2de(Indexbin);
CellSeq = cell(OligoLen,1);
%%
indicesCell = cell(OligoLen, 1);
parfor a = 1:OligoLen
    indicesCell{a} = find(Indexdec == a);
end
for a = 1:OligoLen
    if ~isempty(indicesCell{a})
        CellSeq{a} = DNAmatrix(indicesCell{a}, :);
    else
        CellSeq{a} = []; 
    end
end