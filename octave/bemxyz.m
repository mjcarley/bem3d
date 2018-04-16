function [x,y,z,i,nx,ny,nz]=bemxyz(file)

  ##  [X,Y,Z,I]=bemxyz(FILE): Load a BEM3D geometry dump from FILE.
  ##
  ##   X,Y,Z: node coordinates
  ##   I: mesh node indices
  ##
  ##  To generate the geometry dump from a FILE.BEM
  ##
  ##  bem3d-dump < FILE.BEM > FILE.XYZ

fid = fopen(file, 'r') ;

if ( fid == -1 ) 
  x = y = z = i = nx = ny = nz = [] ;
  warning(["Cannot open file " file]) ;
  return ;
end

dat = fscanf(fid, "%f") ;

fclose(fid) ;

dat = reshape(dat, 7, length(dat)/7)' ;
i = dat(:,1) ;
x = dat(:,2) ; y = dat(:,3) ; z = dat(:,4) ;
[i,j]=sort(i) ;
x = x(j) ; y = y(j) ; z = z(j) ;
nx = dat(j,5) ; ny = dat(j,6) ; nz = dat(j,7) ;

