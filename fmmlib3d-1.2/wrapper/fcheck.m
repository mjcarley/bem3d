k = 1 + j*0.3 ;

load field.dat

x = field(:,1:3) ; 
q = field(:,4) + j*field(:,5) ; 
p = field(:,6) + j*field(:,7) ;

np = length(p) ;

nf = 1024 ;
pp = zeros(nf, 1) ;

for i=1:nf
  ii = [1:i-1 i+1:np] ;

  R = sqrt((x(i,1)-x(ii,1)).^2 + (x(i,2)-x(ii,2)).^2 + (x(i,3)-x(ii,3)).^2) ;

  pp(i) = sum(exp(j*k*R)./R.*q(ii)) ;

endfor
