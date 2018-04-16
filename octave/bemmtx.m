function [A,B]=bemmtx(file)

fid = fopen(file, "r") ;

n = fscanf(fid, "%d", 2) ;
N = n(1) ;

A = B = zeros(N, N) ;

if ( n(2) == 1 ) 
  for i=1:N
    dat = fscanf(fid, "%f", N+1) ; 
    A(dat(1)+1,:) = dat(2:end) ;
    dat = fscanf(fid, "%f", N+1) ; 
    B(dat(1)+1,:) = dat(2:end) ;
  end
else
  for i=1:N
    dat = fscanf(fid, "%f", N+1) ; 
    A(dat(1)+1,:) = dat(2:2:end)+j*dat(3:2:end) ;
    dat = fscanf(fid, "%f", N+1) ; 
    B(dat(1)+1,:) = dat(2:2:end)+j*dat(3:2:end) ;
  end
end

fclose(fid) ;

