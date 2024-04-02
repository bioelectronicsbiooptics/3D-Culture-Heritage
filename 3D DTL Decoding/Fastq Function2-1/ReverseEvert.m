function [teSeq, RevSeqLog] = ReverseEvert(seq,Fsec,Rev,Fw)

Fadd = Fw(1);
if strcmp(Fadd,'A') % strcmp 함수, str 문자열의 일치 비교
    Frp = 'T';
elseif strcmp(Fadd,'T')
    Frp = 'A';
elseif strcmp(Fadd,'C')
    Frp = 'G';
elseif strcmp(Fadd,'G')
    Frp = 'C';
end
%seqrcomplement로 바꿀 수 있는 함수

CompP = seqrcomplement(Rev); % seqrcomplement 함수, revese_complement 시퀀스 출력

S = cell2mat(Fsec); % cell2mat 함수, cell 형식을 matrix 형식으로 변환
F20 = S(:,1:20);
F21 = S(:,2:21);
RevSeqLog = zeros(length(seq),1);
parfor a = 1 : length(seq)
    if sum(F20(a,:) == CompP) >= 18 || sum(F21(a,:) == CompP) >= 18
        RevSeqLog(a) = 1;
    end
end

eSeq = seq(logical(RevSeqLog));
teSeq = cell(1,length(eSeq));
for a = 1 : length(eSeq)
    if sum(eSeq{a}(1:20) == CompP) >= 15
        while length(eSeq{a}) < 150
            eSeq{a} = [eSeq{a} Frp];
        end
        teSeq{a} = seqrcomplement(eSeq{a});
    elseif sum(eSeq{a}(2:21) == CompP) >= 15
        while length(eSeq{a}) < 151
            eSeq{a} = [eSeq{a} Frp];
        end
        teSeq{a} = seqrcomplement(eSeq{a}(2:151));
    end
end