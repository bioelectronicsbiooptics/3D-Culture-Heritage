function matrix = reRSdec(XORM,matrix,k,m,n,Dx,matsz)
d = {'00','01','10','11'};
for a = 1 : length(XORM)
    reRSdata = XORM{a}(:,21:130);
    if ~isempty(reRSdata)
        reSeqtoBin = zeros(size(reRSdata,1),size(reRSdata,2)*2);
        for r = 1 : size(reRSdata,1)
            Nnum = strfind(reRSdata(r,:),'N');
            if  numel(Nnum) >  0
                for i = 1 : numel(Nnum)
                    reRSdata(r,Nnum(i)) = 'A';
                end
            end
            [~,x] = ismember(reRSdata(r,:),'ATGC');
            chr_out = cell2mat(d(x));
            reSeqtoBin(r,:) = chr_out-'0';
        end
        re_payload = reSeqtoBin(:,1:k);
        re_bin_redun = reSeqtoBin(:,k+1:end);
        re_dec_redun = zeros(size(reSeqtoBin,1),(size(reSeqtoBin,2)-k)/m);
        for z = m : m : size(reSeqtoBin,2)-k
            if z == m
                re_dec_redun(:,z/m) = bi2de(re_bin_redun(:,1:m));
            else
                re_dec_redun(:,z/m) = bi2de(re_bin_redun(:,z-m+1:z));
            end
        end
        
        re_rscode = gf([re_payload re_dec_redun],m);
        [re_rxcode, re_cnumerr] = rsdec(re_rscode,n,k);
        re_rxdata = double(re_rxcode.x);
        if ~isempty(find(re_cnumerr==0, 1))
            re_RSdata = mode(re_rxdata(re_cnumerr==0,1:matsz),1);
            matrix(Dx(a),:) = re_RSdata;
        end
    end
end