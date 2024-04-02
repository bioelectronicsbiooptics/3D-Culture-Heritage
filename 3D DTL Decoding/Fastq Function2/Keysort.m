function [oligolen,idxlen,ncp,nbp,ncc,nbc] = Keysort(keymap,KD,k,filename)

strkey = num2str(keymap);
if size(strkey,2) == 5
    ncolposz = str2double(strkey(1));
    nbitposz = str2double(strkey(2));
    ncolcosz = str2double(strkey(3:4));
    nbitcosz = str2double(strkey(5));
    keysz = ncolposz + nbitposz + ncolcosz + nbitcosz;
elseif size(strkey,2) == 6
    ncolposz = str2double(strkey(1:2));
    nbitposz = str2double(strkey(3));
    ncolcosz = str2double(strkey(4:5));
    nbitcosz = str2double(strkey(6));
    keysz = ncolposz + nbitposz + ncolcosz + nbitcosz;
end
d = {'00','01','10','11'};
[~,x] = ismember(KD,'ATGC');
chr_out = cell2mat(d(x));
Keybin = chr_out-'0';

if keymap > 120000
    binkloc = Keybin(length(de2bi(keymap))+9 : length(de2bi(keymap))+9+keysz);
else
    binkloc = Keybin(length(de2bi(keymap))+10:length(de2bi(keymap))+10+keysz);
end
ncp = bi2de(binkloc(1:ncolposz));
nbp = bi2de(binkloc(ncolposz+1:ncolposz+nbitposz));
ncc = bi2de(binkloc(ncolposz+nbitposz+1:ncolposz+nbitposz+ncolcosz));
nbc = bi2de(binkloc(ncolposz+nbitposz+ncolcosz+1:end));

nfile = 3*ncp*nbp+3*ncc*nbc;
tidx = 139;

while tidx + length(de2bi(ceil((nfile)/tidx*1.5))) < k
    tidx = tidx + 1;
end
idxlen = k-tidx;
keyidx = bi2de(Keybin(end-idxlen+1:end));
oligolen = keyidx;
if strcmp(filename,'pisa') || strcmp(filename,'seokga') || strcmp(filename,'liberty')
    oligolen = keyidx+1;
end