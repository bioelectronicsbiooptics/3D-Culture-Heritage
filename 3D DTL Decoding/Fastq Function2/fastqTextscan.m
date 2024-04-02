function seq = fastqTextscan(filename)

for i = 1 : 2
    fid = fopen(filename+"_"+i+".fastq");
    warning('off','all')
    data = textscan(fid,'%s','Delimiter','\n','CommentStyle',{'+','@'});
    warning('on','all')
    if i == 1
        seq_1 = data{1}(2:2:end);
    else
        seq_2 = data{1}(2:2:end);
    end
    fclose(fid);
end
seq = [seq_1;seq_2];