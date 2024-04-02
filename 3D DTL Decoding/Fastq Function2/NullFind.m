function [Ndata,Xdata,dNc,xNc] = NullFind(seq,Nc,Nxc,s,Fp)
% Nc = Null_count; Nxc =Null_xor_count; Fp = F_primers;

if rem(s,2) ~= 0 && rem(Nc,2) ~= 0 && rem(Nxc,2) ~= 0
    dNr = (s/2-Nc/2+1:floor(s/2))+20;
    xNr = (s/2-Nxc/2+1:floor(s/2))+20;
elseif (rem(s,2) ~= 0 && rem(Nc,2) == 0 && rem(Nxc,2) == 0)...
    || (rem(s,2) == 0 && rem(Nc,2) == 0 && rem(Nxc,2) == 0)
    dNr = (floor(s/2)-Nc/2+1:floor(s/2))+20;
    xNr = (floor(s/2)-Nxc/2+1:floor(s/2))+20;
else
    dNr = (floor(s/2)-Nc/2+1:floor(s/2))+20;
    xNr = (floor(s/2)-ceil(Nxc/2)+1:floor(s/2))+20;
end
dpA = repmat('A',1,length(dNr));
xpA = repmat('A',1,length(xNr));

dNo = zeros(length(seq),1);
xNo = zeros(length(seq),1);

%
for a = 1 : length(seq)
    if length(seq{a}) >= s/2+20
        if sum(seq{a}(dNr) == dpA) == length(dpA)
            if sum(seq{a}(1:20) == Fp) >= 15
                dNo(a) = 1;
            end
        end

        if sum(seq{a}(xNr) == xpA) == length(xpA)
            if sum(seq{a}(1:20) == Fp) >= 15
                xNo(a) = 1;
            end
        end
    end
end
%%

dNs = seq(logical(dNo))';
xNs = seq(logical(xNo))';
if gt(numel(xNr),numel(dNr))
    for a = 1 : length(dNs)
        if sum(dNs{a}(xNr) == xpA) == length(xpA)
            dNs{a} = {};
        end
    end
    cellidx = ~cellfun('isempty',dNs);
    dNs = dNs(cellidx);
else
    for a = 1 : length(xNs)
        if sum(xNs{a}(dNr) == dpA) == length(dpA)
            xNs{a} = {};
        end
    end
    cellidx = ~cellfun('isempty',xNs);
    xNs = xNs(cellidx);
end

dNc = cell(length(sum(dNo)),1);
for a = 1 : length(dNs)
    dNc{a} = dNs{a}(1:min(dNr)-1);
end

xNc = cell(length(sum(xNo)),1);
for a = 1 : length(xNs)
    xNc{a} = xNs{a}(1:min(xNr)-1);
end

Ndata = mode(cell2mat(dNc'),1);
Xdata = mode(cell2mat(xNc'),1);