// approximation of a flat plate for scattering calculations
// 
// Parameters:
//    thk: plate thickness;
//    thte: trailing edge angle (12 degrees);
//    len: plate length, not including rounded nose;
//    span: span;
//    msp: span of "metamaterial" section;
//    mle: leading edge of metamaterial section, relative to trailing edge;
//    mlen: length of metamaterial section

thk = 0.05 ;
span = 0.4 ;
len = 0.3 ;
msp = 0.2 ;
mle = 0.2 ;
mlen = 0.1 ;
thte = 12*Pi/180 ;

// Geometric parameters
//

NACA4_th = 12 ; // 12% NACA section for leading edge
NACA4_cut = 0.299827878070145 ;    // point of maximum thickness

Include "naca.geo" ;

NACA4_ch = thk/(NACA4_th/100) ;      // aerofoil chord
NACA4_le_x = 0 ;      // leading edge coordinates
NACA4_le_y = 0 ;      // leading edge coordinates
NACA4_le_z = 0 ;      // leading edge coordinates
NACA4_len_te = 0.009 ; // length scale (trailing edge)
NACA4_len_mp = 0.009 ;  // length scale (mid chord)
NACA4_len_le = 0.008 ;  // length scale (leading edge)
NACA4_nspl = 5   ;    // number of splines on section
NACA4_pps = 5 ;       // number of points per spline
NACA4_cut = 0.299827878070145 ;    // point at which to cut the section
nspl_le = 10 ;         // number of splines for leading edge
nspl_te = 4 ;         // number of splines for trailing edge

// points on section to be used in joining later
p0[] = {} ; p1[] = {} ; p2[] = {} ; p3[] = {} ; p4[] = {} ; p5[] = {} ;

// leading edge points and splines
lpoints[] = {} ; lsplines[] = {} ;
// spanwise and chordwise splines
ssplines[] = {} ; csplines[] = {} ;

// spanwise coordinates of sections
zch[] = {-span/2, -msp/2, msp/2, span/2} ;
NACA4_nspl = nspl_le ; NACA4_le_x = -NACA4_cut*NACA4_ch ;

//NACA4_le_z = -span/2 ; 
For i In {0:3}
  NACA4_le_z = zch[i] ;
  Call NACA4_leading_edge ;
  lsplines[] = {lsplines[], NACA4_Splines[]} ;
  lpoints[] = {lpoints[], NACA4_Points[]} ;
  lle0 = NACA4_Points[0] ; lle1 = NACA4_Points[NACA4_npts-1] ;

  p0[] = {p0[], lle1} ;

  p = newp ; p1[] = {p1[], p} ;
  Point(p) = {len-thk/Tan(thte), -thk/2, NACA4_le_z, NACA4_len_le} ;

  p = newp ; p2[] = {p2[], p} ;
  Point(p) = {len, thk/2, NACA4_le_z, NACA4_len_le} ;

  p = newp ; p3[] = {p3[], p} ;
  Point(p) = {len-mle+mlen, thk/2, NACA4_le_z, NACA4_len_le} ;

  p = newp ; p4[] = {p4[], p} ;
  Point(p) = {len-mle, thk/2, NACA4_le_z, NACA4_len_le} ;

  p5[] = {p5[], lle0} ;
EndFor

// join the sections
For i In {0:2}
  l = newl ; Line(l) = {p0[i], p0[i+1]} ; ssplines[6*i+0] = l ;
  l = newl ; Line(l) = {p1[i], p1[i+1]} ; ssplines[6*i+1] = l ;
  l = newl ; Line(l) = {p2[i], p2[i+1]} ; ssplines[6*i+2] = l ;
  l = newl ; Line(l) = {p3[i], p3[i+1]} ; ssplines[6*i+3] = l ;
  l = newl ; Line(l) = {p4[i], p4[i+1]} ; ssplines[6*i+4] = l ;
  l = newl ; Line(l) = {p5[i], p5[i+1]} ; ssplines[6*i+5] = l ;
EndFor

For i In {0:3}
  l = newl ; Line(l) = {p0[i], p1[i]} ; csplines[5*i+0] = l ;
  l = newl ; Line(l) = {p1[i], p2[i]} ; csplines[5*i+1] = l ;
  l = newl ; Line(l) = {p2[i], p3[i]} ; csplines[5*i+2] = l ;
  l = newl ; Line(l) = {p3[i], p4[i]} ; csplines[5*i+3] = l ;
  l = newl ; Line(l) = {p4[i], p5[i]} ; csplines[5*i+4] = l ;
EndFor

incp = NACA4_nspl*(NACA4_pps - 1)+1 ;

For j In {0:2}
  l0 = ssplines[6*j+5] ;
  For i In {0:(NACA4_nspl-2)}
    l1 = newl ; 
    Line(l1) = {lpoints[j*incp + (NACA4_pps-1)*(i+1)],
      lpoints[(j+1)*incp + (NACA4_pps-1)*(i+1)]} ;
    ll = newll ;
    Printf("Making patch (%g,%g,%g,%g)", 
      l0, lsplines[(j+1)*NACA4_nspl+i], -l1, -lsplines[j*NACA4_nspl+i]) ;    

    Line Loop(ll) = {l0, lsplines[(j+1)*NACA4_nspl+i], -l1, 
      -lsplines[j*NACA4_nspl+i]} ;

    s = news ; 
    Ruled Surface(s) = {-ll} ;
    l0 = l1 ;
  EndFor
  l1 = ssplines[6*j+0] ;
  ll = newll ;
  Printf("Making patch (%g,%g,%g,%g)", 
    l0, lsplines[(j+1)*NACA4_nspl+i], -l1, -lsplines[j*NACA4_nspl+i]) ;    

  Line Loop(ll) = {l0, lsplines[(j+1)*NACA4_nspl+i], -l1, 
    -lsplines[j*NACA4_nspl+i]} ;
  
  s = news ; 
  Ruled Surface(s) = {-ll} ;

EndFor

For j In {0:2}
  For i In {0:4}
    ll = newll ;
    Printf("Making patch (%g,%g,%g,%g)", ssplines[6*j+i], csplines[5*(j+1)+0],
      -ssplines[6*j+i+1], -csplines[5*j+0]) ;
  
		If (j==1 && i == 3)

		Else
			Line Loop(ll) = {ssplines[6*j+i], csplines[5*(j+1)+i], -ssplines[6*j+i+1],
			  -csplines[5*j+i]} ;
			s = news ;
			Ruled Surface(s) = {-ll}; 
		EndIf
  EndFor
EndFor

// tip sections
sgn = -1 ;
For j In {0:3:3}
  l = newl ; Line(l) = {p0[j],p5[j]} ;
  ll = newll ; 
  Line Loop(ll) = {l,lsplines[{j*NACA4_nspl:(j*NACA4_nspl+NACA4_nspl-1)}]} ;
  s = news ; Plane Surface(s) = {sgn*ll} ;

  ll = newll ;
  Line Loop(ll) = {-l,csplines[5*j+0],csplines[5*j+1],csplines[5*j+2],
    csplines[5*j+3],csplines[5*j+4]} ;
  s = news ; Plane Surface(s) = {sgn*ll} ;
  sgn = -sgn ;
EndFor