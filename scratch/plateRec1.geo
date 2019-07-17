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
lp[] = {} ;

// leading edge points and splines
lpoints[] = {} ; lsplines[] = {} ;
// spanwise and chordwise splines
ssplines[] = {} ; csplines[] = {} ;

// spanwise coordinates of sections
zch[] = {-span/2, -msp/2, msp/2, span/2} ;
NACA4_nspl = nspl_le ; NACA4_le_x = -NACA4_cut*NACA4_ch ;

//NACA4_le_z = -span/2 ; 

For i In {1:2}
  NACA4_le_z = zch[i] ;

  p = newp ; p3[] = {p3[], p} ;
  Point(p) = {len-mle+mlen, thk/2, NACA4_le_z, NACA4_len_le} ;

  p = newp ; p4[] = {p4[], p} ;
  Point(p) = {len-mle, thk/2, NACA4_le_z, NACA4_len_le} ;

EndFor

//+
l =newl; lp[] = {lp[], l} ; Line(l) = {p4[0], p3[0]} ;
//+
l =newl; lp[] = {lp[], l} ; Line(l) = {p3[0], p3[1]} ;
//+
l =newl; lp[] = {lp[], l} ; Line(l) = {p3[1], p4[1]} ;
//+
l =newl; lp[] = {lp[], l} ; Line(l) = {p4[1], p4[0]} ;
//+
Line Loop(1) = {lp[0], lp[1], lp[2], lp[3]};
//+
Plane Surface(1) = {-1};
