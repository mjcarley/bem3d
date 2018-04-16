function [f, i]=bemdat(file)

  ##  [F,I]=bemdat(FILE): Load a BEM3D mesh data block from FILE.
  ##
  ##   F: block of data, each row corresponding to a mesh node
  ##   I: sorted mesh node indices. 

fid = fopen(file, 'r') ;

if ( fid == -1 ) 
  i = f = [] ;
  warning(["Cannot open file " file]) ;
  return ;
end

n = fscanf(fid, '%d', 2) ;
w = n(2) ; n = n(1) ;
s = fscanf(fid, '%s', 1) ;

dat = fscanf(fid, '%f') ;

i = dat(1:(w+1):end) ;
dat(1:(w+1):end) = [] ;
f = reshape(dat,w,n)' ;

fclose(fid) ;

[i,j] = sort(i) ;
f = f(j,:) ;

