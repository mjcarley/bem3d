function p=isphere(k, a, r, th, bt, M)

  ## P = ISPHERE(K, A, R, THETA, BETA)
  ##
  ## Compute the total acoustic field at spherical coordinates (R,
  ## THETA) of wavenumber K around a sphere of radius A and constant
  ## admittance BETA, irradiated by a plane wave propagating in the
  ## positive x direction, using the formulae from Morse and Ingard,
  ## Theoretical Acoustics, p425. Note that Morse and Ingard define the
  ## admittance such that dp/dn = -j*k*BETA*p on the surface

if ( nargin < 5 ) bt = 0 ; endif
if ( nargin < 6 ) M = 10 ; endif

C = cos(th) ;
ka = k*a ;
kr = k*r ;

Pm = ones(size(th)) ; 
Pmp1 = C ;

p = exp(j*k*r.*C) ;

## using the Bessel/Hankel recursions keeps the function
## self-contained

##jm = sin(ka)/ka ; jmp1 = sin(ka)/ka^2 - cos(ka)/ka ;
##hm = exp(j*ka)/(j*ka) ; hmp1 = exp(j*ka)/ka*(1/(j*ka)-1) ;

## (but it's unstable so I went back to using the built-in versions)

##hmr = exp(j*kr)./(j*kr) ; hmp1r = exp(j*kr)./kr.*(1./(j*kr)-1) ;

for m=0:M
  hm = sqrt(pi/2./ka).*besselh(m+0.5, 1, ka) ;

  hp1 = sqrt(pi/2./ka).*besselh(m+1+0.5, 1, ka) ;
  dhm = -hp1 + m*hm./ka ;

  jm = real(hm) ; djm = real(dhm) ;
  hmr = sqrt(pi/2./kr).*besselh(m+0.5, 1, kr) ;

  tt = (djm + j*bt*jm)./(dhm + j*bt*hm) ;
  p -= (2*m+1)*j^m*Pm.*tt.*hmr ;

  tt = Pm ; Pm = Pmp1 ;
  Pmp1 = (2*m+3)*C.*Pm - (m+1)*tt ;
  Pmp1 /= m+2 ;

endfor
