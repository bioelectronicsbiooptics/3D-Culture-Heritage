function [tc,tp] = tridecode(ECCp,ncp,nbp,ncc,nbc,filename)

if rem(ncp*3*nbp+ncc*3*nbc,2) ~=0
    ECCp(end) = [];
end

pr = ECCp(1:ncp*nbp*3);
cr = ECCp(ncp*nbp*3+1:end);

pb = reshape(pr',[nbp ncp*3'])';
pd = bi2de(double(pb))/(10^16);

tp = reshape(pd,[3 ncp])';
if strcmp(filename,'seokga')||strcmp(filename,'liberty')||strcmp(filename,'pisa')...
        ||strcmp(filename,'venus')||strcmp(filename,'goldengate')||strcmp(filename,'tripitaka')
    tp = reshape(pd,[ncp 3]);
end

cb = reshape(cr',[nbc ncc*3])';
cd = bi2de(double(cb));
tc = reshape(cd,[3 ncc])';
% ERROR
tc(find(tc==0)) = 1;
tc(find(tc>length(tp))) = 1;