function p=sphscat(a,r,th,k,b)

  ## P=SPHSCAT(A,R,TH,K): Acoustic scattering by a sphere
  ##
  ## Calculate the total acoustic field around a sphere subject to a
  ## plane wave. 
  ## 
  ## A:  sphere radius
  ## R:  field position radius (a vector)
  ## TH: field position azimuthal angle, measured from propagation
  ##     incident wave
  ## K:  wavenumber
  ##
  ## P=SPHSCAT(A,R,TH,K,B) calculates the field for a sphere with surface
  ##   admittance B

if ( nargin < 5 ) b = 0 ; end

if ( any(size(r) ~= size(th)) )
  if ( size(r,1) == 1 & size(r,2) == 1 )
    r = r*ones(size(th)) ;
  end
  if ( size(th,1) == 1 & size(th,2) == 1 )
    th = th*ones(size(r)) ;
  end
  if ( any(size(r) ~= size(th)) )
    error("Coordinate matrices R and THETA must have conforming sizes.") ;
  end
end

r = r(:)' ; th = th(:)' ;

p = zeros(size(r)) ;
ii = find(r >= a) ;

for m=0:40
  [jm,djm] = sbessel1(m,k*a) ; [nm,dnm] = sbessel2(m,k*a) ; 
  hm = jm + j*nm ; dhm = djm + j*dnm ;
  Rm = (conj(dhm) + j*b*conj(hm))/(dhm + j*b*hm) ;
  Pm = legendre(m,cos(th(ii))) ; Pm = Pm(1,:) ;
  jm = sbessel1(m,k*r(ii)) ; nm = sbessel2(m,k*r(ii)) ; 
  hm = jm + j*nm ;
  p(ii) += (2*m+1)*j^m*Pm.*(jm-0.5*(1+Rm)*hm) ;
end
