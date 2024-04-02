function [teSeq, RevSeqLog] = RecCompOrder(seq,Fsec,Rev,Fw)
% Fsec = ForwardSection;
% Rev = R_primers;
% Fw = F_primers;

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

% %% Option 1 (T-21) 
% parfor a = 1:length(eSeq)
%     seqLength = length(eSeq{a});
% 
%     % 1:20 문자열 처리
%     if sum(eSeq{a}(1:20) == CompP) >= 15 && seqLength >= 150
%         Reverse = eSeq{a}(1:150);
%         teSeq{a} = seqrcomplement(Reverse);
%     end
% 
%     % 2:21 문자열 처리
%     if sum(eSeq{a}(2:21) == CompP) >= 15 && seqLength >= 151
%         Reverse = eSeq{a}(2:151);
%         teSeq{a} = seqrcomplement(Reverse);
%     end
% end
% %% Option 2 (S-21) 
% parfor a = 1:length(eSeq)
%     seqLength = length(eSeq{a});
%     extendLength = max(150 - seqLength, 0); % 150보다 작을 경우 확장해야 할 길이 계산
% 
%     if sum(eSeq{a}(1:20) == CompP) >= 15 && seqLength < 150
%         eSeq{a} = [eSeq{a} repmat(Frp, 1, extendLength)]; % 한 번에 확장
%         teSeq{a} = seqrcomplement(eSeq{a}(1:150));
%     elseif sum(eSeq{a}(2:21) == CompP) >= 15 && seqLength < 151
%         extendLength = max(151 - seqLength, 0); % 151보다 작을 경우 확장해야 할 길이 계산
%         eSeq{a} = [eSeq{a} repmat(Frp, 1, extendLength)]; % 한 번에 확장
%         teSeq{a} = seqrcomplement(eSeq{a}(2:151));
%     end
% end
%% 39?
parfor a = 1:length(eSeq)
    seqLength = length(eSeq{a});
    extendLength = 150 - seqLength; % 150보다 작을 경우 확장해야 할 길이 계산

    if sum(eSeq{a}(1:20) == CompP) >= 15
        eSeq{a} = [eSeq{a} repmat(Frp, 1, extendLength)]; % 한 번에 확장
        teSeq{a} = seqrcomplement(eSeq{a});
    elseif sum(eSeq{a}(2:21) == CompP) >= 15
        extendLength = 151 - seqLength; % 151보다 작을 경우 확장해야 할 길이 계산
        eSeq{a} = [eSeq{a} repmat(Frp, 1, extendLength)]; % 한 번에 확장
        teSeq{a} = seqrcomplement(eSeq{a}(2:151));
    end
end
%%
% parfor a = 1 : length(eSeq) % 608sec (normal) 113sec (parfor)
%     if sum(eSeq{a}(1:20) == CompP) >= 15
%         while length(eSeq{a}) < 150
%             eSeq{a} = [eSeq{a} Frp];
%         end
%         teSeq{a} = seqrcomplement(eSeq{a});
%     elseif sum(eSeq{a}(2:21) == CompP) >= 15
%         while length(eSeq{a}) < 151
%             eSeq{a} = [eSeq{a} Frp];
%         end
%         teSeq{a} = seqrcomplement(eSeq{a}(2:151));
%     end
% end
