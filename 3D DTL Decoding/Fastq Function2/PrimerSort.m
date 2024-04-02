function Pseq = PrimerSort(seq,F,R)

Pseq = cell(length(seq),1);
for a = 1 : length(seq)
    temp = seq{a};
     if length(temp) >= 150
        if sum(temp(1:20) == F) >=18
            if sum(temp(131:150) == R) >= 18
                Pseq{a} = temp(1:150);
            end
        end
     end
    if length(temp) >= 151
        if sum(temp(2:21) == F) >=18
            if sum(temp(132:151) == R) >= 18
                Pseq{a} = temp(2:151);
            end
        end
    end
end
Pseq = Pseq(~cellfun('isempty',Pseq));
