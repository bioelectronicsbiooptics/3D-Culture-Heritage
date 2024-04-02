function [key4dec, inf_key] = genkey(stl_points, stl_connects, bin_points, bin_connects)

ncolpo = length(stl_points); 
ncolco = length(stl_connects);
nbitpo = size(bin_points,2); 
nbitco = size(bin_connects,2);

inf_stl.cp = (de2bi(ncolpo));
inf_stl.cc = (de2bi(ncolco));
inf_stl.bp = (de2bi(nbitpo));
inf_stl.bc = (de2bi(nbitco));

inf_key = [inf_stl.cp inf_stl.bp inf_stl.cc inf_stl.bc];

key4dec = length(inf_stl.cp)*10^4+length(inf_stl.bp)*10^3+length(inf_stl.cc)*10+length(inf_stl.bc);