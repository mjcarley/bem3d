function [h,dh]=shankel(n, x)

h = sqrt(pi/2./x).*besselh(n+0.5, 1, x) ;

if ( nargout < 2 ) return ; endif

hp1 = sqrt(pi/2./x).*besselh(n+1+0.5, 1, x) ;
dh = -hp1 + n*h./x ;
