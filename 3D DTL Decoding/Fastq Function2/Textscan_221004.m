%%
addpath Tripitaka Tower Jikji V&G
filename = input("filename?: ",'s');
tic
for i = 1 : 2
    fid = fopen(filename+"_"+i+".fastq");
    data = textscan(fid,'%s','Delimiter','\n','CommentStyle',{'+','@'});
    if i == 1
        seq_1 = data{1}(2:2:end);
    else
        seq_2 = data{1}(2:2:end);
    end
    fclose(fid);
end
seq = [seq_1;seq_2];
toc