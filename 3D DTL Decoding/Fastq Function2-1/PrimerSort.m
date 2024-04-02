function Pseq = PrimerSort(seq,F,R)

Pseq = cell(length(seq),1);
parfor a = 1 : length(seq)
    temp = seq{a};

    if length(temp) >= 150 && sum(temp(1:20) == F) == 20 && sum(temp(131:150) == R) == 20
                Pseq{a} = temp(1:150);
    end
    if length(temp) >= 151 && sum(temp(2:21) == F) == 20 && sum(temp(132:151) == R) == 20
                Pseq{a} = temp(2:151);
    % % elseif length(temp) >= 151
        if sum(temp(2:21) == F) >=18
            if sum(temp(132:151) == R) >= 18
                Pseq{a} = temp(2:151);
            end
        end
    end

    W = seq{a};
    if length(W) >= 150
        if sum(W(1:20) == F) >=18
            if sum(W(131:150) == R) >= 18
                Pseq{a} = W(1:150);
            end
        end
    end
    if length(W) >= 151
        if sum(W(2:21) == F) >=18
            if sum(W(132:151) == R) >=18
                Pseq{a} = W(2:151);
            end
        end
    % elseif length(temp) >= 151
        if sum(temp(2:21) == F) >=18
            if sum(temp(132:151) == R) >= 18
                Pseq{a} = temp(2:151);
            end
        end
    end

end
Pseq = Pseq(~cellfun('isempty',Pseq)); % cellfun 함수, cell 형식 변수의 각 cell 마다 함수 적용/ isempty 함수, 비어 있는지 확인하여 논리값 출력
% size(Pseq,1)