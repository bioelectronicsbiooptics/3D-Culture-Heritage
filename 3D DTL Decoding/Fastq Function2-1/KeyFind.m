function [Key,KmatDNA] = KeyFind(PrimerSortingSeq)

DNAmatrix = cell2mat(PrimerSortingSeq);
payload = DNAmatrix(:,21:98);
KeyM = zeros(length(payload),1);
parfor a = 1 : size(payload,1)
    if numel(strfind(payload(a,:),'AAAAA')) >= 1
        if numel(strfind(payload(a,:),'AAAA')) - ...
                numel(strfind(payload(a,:),'AAAAA')) + ...
                numel(strfind(payload(a,:),'GAAAT')) >= 3 %?
            KeyM(a) = 1;
        end
    end
end
pK = PrimerSortingSeq(find(KeyM),:);
mK = cell2mat(pK);

Ks = cell(size(mK,1),1);
for i = 1 : size(mK,1)
    Ks{i} = mK(i,21:98);
end

mDNA = [repmat('TGCA',1,27) 'TG'];
KeyD = cell(length(Ks),1);
for a = 1 : length(Ks)
    temp = Ks{a};
    for b = 73 : -1 : 58 % index 앞자리
        keyi = length(temp(b+1:end));%index 길이 = poly a 길이
        if strcmp(temp(b-4:b),repmat('A',1,5)) % poly a 5개가 있는지 확인하는 것
            if strcmp(temp(b-keyi-10:b-keyi-1),mDNA(b-keyi-10:b-keyi-1)) %null을 찾는 것
                KeyD{a} = temp;
            end
        end
    end
end
KeyD = KeyD(~cellfun('isempty',KeyD));

KmatDNA = mode(cell2mat(KeyD),1); %최빈값 각 자리로 하는지 통으로 하는지 확인해보기.
d = {'00','01','10','11'};
for b = 1 : size(KmatDNA,2)-3
    if strcmp(KmatDNA(b:b+3),'AAAA')
        KeyMaptemp = KmatDNA(1:b-1);
        if (strfind(KeyMaptemp,'N'))>0
            KeyMaptemp = repmat('A',1,length(KeyMaptemp));%?
        end
        [~,x] = ismember(KeyMaptemp,'ATGC');
        chr_out = cell2mat(d(x));
        tempKeyMap = chr_out-'0';
        if ~isempty(tempKeyMap)
            Key = bi2de(tempKeyMap);
        end
        break
    end
end