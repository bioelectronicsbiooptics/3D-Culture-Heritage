function [c,x,resh_data,index,Maplocation] = cCal(column_data,xor_data,k,filename)

heritg = categorical({'seokga' 'liberty' 'pisa'}); 
fileloc = (heritg == filename);
Lat = []; Lng = [];
if fileloc(1) == 1
    Lat = 35.789997826271744;
    Lng = 129.33202224714583;
elseif fileloc(2) == 1
    Lat = 40.68926957444817;
    Lng = -74.04461479091522;
    Lng = Lng + 110;
elseif fileloc(3) == 1
    Lat = 43.72295857957317;
    Lng = 10.396633498022156;    
    Lng = Lng + 20;
end

if isempty(Lat) == 0 
    binLat = de2bi((Lat+150)*10^14);
    binLng = de2bi((Lng+150)*10^14);
end

Null_count = 0;
while rem(numel(column_data),k) ~= 0
    column_data(end+1) = 0;
    Null_count = Null_count + 1;
end

c = Null_count;
rd = reshape(column_data',[k numel(column_data)/k])';

Null_xor_count = 0;
while rem(numel(xor_data),k) ~= 0
    xor_data(end+1) = 0;
    Null_xor_count = Null_xor_count + 1;
end


x = Null_xor_count;
rx = reshape(xor_data',[k numel(xor_data)/k])';

resh_data = [rd; rx];

oligo_num = size(resh_data,1);

index = de2bi(1 : oligo_num);
Maplocation = [];
if isempty(Lat) == 0
    mapidx = de2bi(oligo_num + 1);
    mapNull = [0 0 0 1 1 0 1 1];
    mapdata = [binLat binLng];

    while length(mapdata) < 156 - size(mapidx,2)
        mapdata = [mapdata mapNull];
    end

    while length(mapdata) > 156 - size(mapidx,2)
        mapdata(end) = [];
    end

    mapdata = [mapdata mapidx];

    mk = 156;
    mn = mk + 8;
    mm = length(de2bi(mn));

    gfmap = gf(mapdata,mm); 
    rsmap = rsenc(gfmap,mn,mk); 

    rsmap = rsmap.x; 
    binrsmap = rsmap(1:mk);
    rdcrsmap = rsmap(mk+1:end);
    rdbit = de2bi(rdcrsmap);

    rsmapdata = [binrsmap reshape(rdbit,[1 numel(rdbit)])];

    maploc = logical(rsmapdata); 
    bases = 'ATGC';
    for a = 1 : 2 : length(maploc)   

      DNA_index = 2 * maploc(1,a) + maploc(1,a+1) + 1; 

      DNAmaploc(1,(a+1)/2) = bases(DNA_index);

    end
    Maplocation = [Maplocation DNAmaploc];
end