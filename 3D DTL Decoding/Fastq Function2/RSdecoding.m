function [rx,err] = RSdecoding(bS,k,m,n)

bP = bS(:,1:k);
bR = bS(:,k+1:end);
dR = zeros(size(bS,1),(size(bS,2)-k)/m);
for z = m : m : size(bS,2)-k
    if z == m 
        dR(:,z/m) = bi2de(bR(:,1:m));
    else
        dR(:,z/m) = bi2de(bR(:,z-m+1:z));
    end
end

rscode = gf([bP dR],m);
if length(rscode)<1000
    [rxcode, err] = rsdec(rscode,n,k);
    rxcode = rxcode.x;
    rx = double(rxcode);
else
    mR = gf(zeros(size(rscode,1),size(rscode,2)-(n-k)),m);
    err = zeros(size(rscode,1),1);
    for a = 1000 : 1000 : size(rscode,1)
        [rxcode, cnumerr] = rsdec(rscode(a-999:a,:),n,k);
        mR(a-999:a,:) = rxcode;
        err(a-999:a) = cnumerr;
    end
    [rxcode, cnumerr] = rsdec(rscode(a:end,:),n,k);
    mR(a:end,:) = rxcode;
    err(a:end,:) = cnumerr;
    mR = mR.x;
    rx = double(mR);
end