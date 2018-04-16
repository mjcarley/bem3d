function i=bemedge(file)

##  I=bemedge(FILE): Load a BEM3D edge file from FILE
##
##  I: node indices for a sharp edge on a BEM3D mesh.

fid = fopen(file, "r") ;

if ( fid == -1 )
  error(["Cannot open " file " for reading"]) 
end

ne = fscanf(fid, "%d", 1) ;
fscanf(fid, "%[^\n]c") ;
fscanf(fid, "%c", 1) ;

[i,ni] = fscanf(fid, "%d", 2*ne) ;

if ( ni ~= 2*ne )
  error(["File " file " should have " int2str(2*ne) "indices and only has " ...
	 int2str(ni)]) ;
end

i = [i(1:2:end) i(2:2:end)] ;

fclose(fid) ;