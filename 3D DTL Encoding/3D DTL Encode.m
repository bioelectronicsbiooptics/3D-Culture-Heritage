% STL_Encode
%% STL ENCODING 
%% STL File open

filename = input("file name:",'s'); 
stlid = stlread(filename+".stl"); 

stl_connects = stlid.ConnectivityList; 
bin_connects = de2bi(stl_connects'); 
stl_points = stlid.Points; 
% if  ~regexp(filename, regexptranslate('wildcard','Jikji*')) 
%    if min(stl_points) < 0
        stl_points = stl_points + min(stl_points)*(-1);
%    end
% end
bin_points = de2bi((stl_points+1)*10^16); 
%% Key
[key4dec, inf_key] = genkey(stl_points, stl_connects, bin_points, bin_connects); 

po2col = reshape(bin_points',[numel(bin_points) 1])'; 
co2col = reshape(bin_connects',[numel(bin_connects) 1])';
column_data =[po2col co2col];

if rem(numel(column_data),2) ~= 0 
    column_data(end+1) = 0;
end
%% XOR
data_half = numel(column_data)/2;
xor_A = column_data(1:data_half);
xor_B = column_data(data_half+1:end);
xor_data = xor(xor_A,xor_B); 
%% index
k = 139;
[Null_col,Null_xor,resh_data,index,Maplocation] = cCal(column_data,xor_data,k,filename); % 첨부된 cCal 함수 참조

while size(index,2) + k < 156 
    k = k + 1;
    [Null_col,Null_xor,resh_data,index,Maplocation] = cCal(column_data,xor_data,k,filename);
end
%% RS code
k = k + size(index,2);
n = k + 8 ;
m = length(de2bi(n));

nrsin = [resh_data index];
nrsin(end+1,:) = 0;

key2seq_in = de2bi(key4dec);
Lkeyseq = length(key2seq_in);

nrsin(end,1:Lkeyseq) = key2seq_in;

if rem(Lkeyseq,2) ~= 0
    nrsin(end,Lkeyseq+10:Lkeyseq+9+length(inf_key)) = inf_key;
    key_end = Lkeyseq+9+length(inf_key);
else
    nrsin(end,Lkeyseq+9:Lkeyseq+8+length(inf_key)) = inf_key;
    key_end = Lkeyseq+8+length(inf_key);
end

nrsin(end,end-size(index,2)+1:end) = de2bi(length(index)+1);
nNull = [de2bi(Null_col) de2bi(Null_xor)];
nrsin(end,key_end+9:key_end+8+length(nNull)) = nNull;

gf_en = gf(double(nrsin), m);
code = rsenc(gf_en, n, k);

gf_payload = code(:,1:k);
gf_redundancy = code(:,k+1:n);
bin_payload = gf_payload.x;
dec_redun = gf_redundancy.x';

for i = 1 : size(code,1)
    dr = de2bi(dec_redun(:,i));
    while size(dr,2) < m
        dr(:,end+1) = 0;
    end
    br = reshape(dr',[1 numel(dr)]);
    bin_redun(i,1:m*(n-k)) = br;
end

bin_seq = [bin_payload bin_redun];
%% Encode to DNA
lo_rscode = logical(bin_seq); 
bases = 'ATGC';

for r = 1 : size(bin_seq,1)
    for a = 1 : 2 : size(lo_rscode,2)   

      DNA_index = 2 * lo_rscode(r,a) + lo_rscode(r,a+1) + 1; 

      result_code(r,(a+1)/2) = bases(DNA_index);
    end
end

as_DNA_code = result_code;
for c = 1: size(result_code,2)
    mode = rem(c,4) + 1;
    modulation(1,c) = bases(mode);
end

as_DNA_code(end,(key_end+1)/2+15:end-(size(bin_redun,2)/2+size(index,2)))...
= modulation((key_end+1)/2+15:end-(size(bin_redun,2)/2+size(index,2)));

[F_primers, R_primers] = primerSelect(filename); 
%% Save as DNA file
for d = 1 : 1 : size(as_DNA_code,1)
    addition = [F_primers as_DNA_code(d,:) R_primers];
    DNA_library(d,:) = addition;
end

if isempty(Maplocation) == 0
    DNA_library = [DNA_library; ([F_primers Maplocation R_primers])];
end
writematrix(DNA_library,filename+".txt"); 