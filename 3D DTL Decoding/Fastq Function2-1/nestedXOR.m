function [XORM, Dx] = nestedXOR(CellSeq,xtp,ErrD,RSerrinf,XA,XB,XC,payload,tp,xp,matsz,Null_count)

Depth_cut = ErrD(RSerrinf(ErrD)~=0);
A_part = Depth_cut(Depth_cut<=size(xp,1));
C_part = Depth_cut(Depth_cut>size(tp,1));
B_part = Depth_cut(length(A_part)+1:end-length(C_part));

if ~isempty(B_part)
    if B_part(end) == size(tp,1)
        B_part(end) = [];
    end
end

if ~isempty(C_part)
    if C_part(end) == size(tp,1)+size(xp,1)-1
        C_part(end) = [];
    end
end

Xi = zeros(length(Depth_cut),1);
for a = 1 : length(A_part)
    AtoB = [((A_part(a)-1)*matsz)+1 ((A_part(a)-1)*matsz)+matsz]- ...
        rem(length(XA),matsz);
    AtoC = [(A_part(a)-1)*matsz+1 (A_part(a)-1)*matsz+matsz];
    if sum(XB(AtoB(1):AtoB(2)))==0 || sum(XC(AtoC(1):AtoC(2)))==0
        Xi(a) = A_part(a);
    end
end

for b = 1 : length(B_part)
    BtoA = [((B_part(b))*matsz)+1 ((B_part(b))*matsz)+matsz]- ...
        length(XA)-(matsz-rem(length(XA),matsz));
    BtoC = [((B_part(b))*matsz)+1 ((B_part(b))*matsz)+matsz]- ...
        length(XA)-(matsz-rem(length(XA),matsz));
    if sum(XA(BtoA(1):BtoA(2)))==0 || sum(XC(BtoC(1):BtoC(2)))==0
        Xi(b+length(A_part)) = B_part(b);
    end
end

for c = 1 : length(C_part)
    CtoA = [(C_part(c)-1)*matsz+1 (C_part(c)-1)*matsz+matsz]- ...
        length(XA)*2-Null_count;
    CtoB = [((C_part(c)-1)*matsz)+1 ((C_part(c)-1)*matsz)+matsz]- ...
        length(XA)*2-rem(length(XA),matsz)-Null_count;
    if sum(XB(CtoA(1):CtoA(2)))==0 || sum(XB(CtoB(1):CtoB(2)))==0
        Xi(c+length(A_part)+length(B_part)) = C_part(c);
    end
end

Non_XOR_oligo=Xi(find(Xi));

for i = 1 : length(Non_XOR_oligo)
    if RSerrinf(Non_XOR_oligo(i)) == 0 && Non_XOR_oligo(i)<size(xtp,1)
        xtp(Non_XOR_oligo(i),:) = payload(Non_XOR_oligo(i),:);
    elseif RSerrinf(Non_XOR_oligo(i)) == -1
        
    end
end

Dx = Non_XOR_oligo(find(RSerrinf(Non_XOR_oligo)));
XORM = CellSeq(Dx);