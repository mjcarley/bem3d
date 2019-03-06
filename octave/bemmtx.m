function [A,B]=bemmtx(file)

  ## [A,B]=bemmtx(FILE): read BEM3D's solution matrices from FILE
  ##
  ## Matrices read are such that [A]phi = [B]dphi
  ## Currently only implemented for direct solver (with dense matrices)
  
fid = fopen(file, "r") ;

solver = fscanf(fid, "%[^ ]c") ;
##solver

n = fscanf(fid, "%d", 6) ;
N = n(2) ;
w = n(3) ;

##[N w]

A = B = zeros(N, N) ;

if ( w == 1 ) 
  for i=1:N
    dat = fscanf(fid, "%f", N+1) ; 
    A(dat(1)+1,:) = dat(2:end) ;
    dat = fscanf(fid, "%f", N+1) ; 
    B(dat(1)+1,:) = dat(2:end) ;
  end
else
  for i=1:N
    dat = fscanf(fid, "%f", 2*N+1) ;
    ##    size(dat)
    ##    dat(1)
    A(dat(1)+1,:) = dat(2:2:end)+j*dat(3:2:end) ;
    dat = fscanf(fid, "%f", 2*N+1) ; 
    B(dat(1)+1,:) = dat(2:2:end)+j*dat(3:2:end) ;
  end
end

fclose(fid) ;
