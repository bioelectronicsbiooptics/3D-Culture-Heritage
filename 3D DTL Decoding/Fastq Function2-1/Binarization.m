function binSeq = Binarization(CellSeq)

dataMode = cell(length(CellSeq),1);
for a = 1 : length(CellSeq)
    if ~isempty(CellSeq{a})
        dataMode{a} = mode(CellSeq{a},1); % mode 함수, 최빈값
    end
end
dataMode(cellfun('isempty',dataMode)) = {repmat('A',1,150)};
DNAdata = cell2mat(dataMode);

MatData = DNAdata(:,21:130);
binSeq = zeros(length(MatData),size(MatData,2)*2); % length 잘못된걸수도. size(MatData,1)로
d = {'00','01','10','11'};
for r = 1 : size(MatData,1)
    Nnum = strfind(MatData(r,:),'N');
    if  numel(Nnum) >  0
        for i = 1 : numel(Nnum)
            MatData(r,Nnum(i)) = 'A';
        end
    end
    [~,x] = ismember(MatData(r,:),'ATGC');
    chr_out = cell2mat(d(x));
    binSeq(r,:) = chr_out-'0';
end