k = 3.9 ; a = 1 ;

[x,y,z,i]=bemxyz("plane.xyz");
[i,f]=bemdat("plane.dat") ;

R = sqrt(x.^2+y.^2+z.^2)' ;
th = atan2(y,x)' ;

p = sphscat(a,R,th,k) ;

q = f(:,1) + j*f(:,2) ;

ip = exp(j*k*x') ;
