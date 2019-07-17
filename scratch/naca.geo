// for details of computation methods, see NASA TM 4741, Ladson,
// Brooks, Hill and Sproles, 1996.

// constants from NASA TM4741
NACA4_a0 =  0.2969 ;
NACA4_a1 = -0.1260 ;
NACA4_a2 = -0.3516 ;
NACA4_a3 =  0.2843 ;
NACA4_a4 = -0.1015 ;

Function NACA4

// Generate a symmetric NACA four series aerofoil composed of a set of
// splines. The section is modified to have the trailing edge extended
// slightly to form a true sharp trailing edge, rather than the blunt edge
// of the real section. This means that a section of nominal unit chord is
// actually slightly longer (1.015 rather than 1, in the case of a 0020 
// section).
//
// the following variables must be set on entry:
// NACA4_th     thickness in percent of chord
// NACA4_ch     aerofoil chord
// NACA4_le_x,y,z   leading edge coordinates
// NACA4_len_te length scale (trailing edge)
// NACA4_len_mc length scale (mid chord)
// NACA4_len_le length scale (leading edge)
// NACA4_nspl   number of splines on section
// NACA4_pps    number of points per spline

// The local scale length will be set using a quadratic which interpolates
// trailing edge, midpoint and leading edge scale lengths as a function of
// distance along the chord.

// On exit, the following variables will contain the details of
// a closed NACA 00TH section with 2*NACA4_nspl splines, each of
// NACA4_pps points:
//
// NACA4_Points[] a list of the 2*NACA4_nspl*NACA4_pps points
//                on the section
// NACA4_Splines[] a list of the 2*NACA4_nspl splines
//
// these two lists are oriented so that they start at the trailing edge
// and move over the upper surface and around the lower surface to return
// to the trailing edge

NACA4_npts = NACA4_nspl*(NACA4_pps-1)-1 ;

NACA4_Points[] = {} ;
NACA4_Splines[] = {} ;

// Generate the trailing edge extension
_x = 1 ;
_y = NACA4_a0*Sqrt(_x) + 
    _x*(NACA4_a1 + _x*(NACA4_a2 + _x*(NACA4_a3 + _x*NACA4_a4))) ;
_dy = 0.5*NACA4_a0/Sqrt(_x) + 
    (NACA4_a1 + _x*(2*NACA4_a2 + _x*(3*NACA4_a3 + _x*4*NACA4_a4))) ;

_x = 1 - _y/_dy ; _y = 0 ;
_p = newp ;
Point(_p) = {NACA4_ch*_x+NACA4_le_x,  
	    _y*NACA4_ch*NACA4_th/20+NACA4_le_y, 
	    NACA4_le_z, NACA4_len_te} ;
NACA4_Points[0] = _p ;

For _i In {1:NACA4_npts}
  _th = (_i-1)/(NACA4_npts-1)*2*Pi ;
  _x = 0.5*(1+Cos(_th)) ;
  _y = NACA4_a0*Sqrt(_x) + 
  _x*(NACA4_a1 + _x*(NACA4_a2 + _x*(NACA4_a3 + _x*NACA4_a4))) ;

  If ( _th > Pi ) _y = -_y ; EndIf
  
  _p = newp ;
  L1 = 1-_x ; L2 = _x ;
  NACA4_len = NACA4_len_le*L1*(2*L1-1) + NACA4_len_mp*4*L1*L2 +
  NACA4_len_te*L2*(2*L2-1) ;
  Point(_p) = {NACA4_ch*_x+NACA4_le_x,  
    _y*NACA4_ch*NACA4_th/20+NACA4_le_y, 
    NACA4_le_z, NACA4_len} ;
  NACA4_Points[_i] = _p ;
EndFor

For _i In {0:NACA4_nspl-2}
  _c = newc ;
  Spline(_c) = NACA4_Points[{_i*(NACA4_pps-1):(_i+1)*(NACA4_pps-1)}] ;
  NACA4_Splines[_i] = _c ;
EndFor
_c = newc ; _i = NACA4_nspl - 1 ;

Spline(_c) = {NACA4_Points[{_i*(NACA4_pps-1):(NACA4_npts-1)}],
	     NACA4_Points[0]} ;

NACA4_Splines[_i] = _c ;

NACA4_npts += 1 ;

Return

Function NACA4_leading_edge

// Generate the leading edge of a symmetric NACA four series aerofoil
// composed of a set of splines, starting from a specified percentage
// of chord

//
// the following variables must be set on entry:
// NACA4_th     thickness in percent of chord
// NACA4_ch     aerofoil chord
// NACA4_le_x,y,z   leading edge coordinates
// NACA4_len_te length scale (trailing edge)
// NACA4_len_mc length scale (mid chord)
// NACA4_len_le length scale (leading edge)
// NACA4_nspl   number of splines (best if this is odd)
// NACA4_pps    number of points per spline
// NACA4_cut    point at which to cut the section, expressed as a 
//              fraction of chord

// The local scale length will be set using a quadratic which interpolates
// trailing edge, midpoint and leading edge scale lengths as a function of
// distance along the chord.

// On exit, the following variables will contain the details of
// the leading edge
//
// NACA4_Points[] a list of the 2*NACA4_nspl*NACA4_pps points
//                on the section
// NACA4_Splines[] a list of the 2*NACA4_nspl splines
//
// these two lists are oriented so that they start at the trailing edge
// and move over the upper surface and around the lower surface to return
// to the trailing edge

NACA4_npts = NACA4_nspl*(NACA4_pps-1) + 1 ;

NACA4_Points[] = {} ;
NACA4_Splines[] = {} ;

_th0 = Acos(2*NACA4_cut - 1) ;

For _i In {0:NACA4_npts-1}
  _th = _th0 + _i/(NACA4_npts-1)*(2*Pi - 2*_th0) ;
  _x = 0.5*(1+Cos(_th)) ;
  _y = NACA4_a0*Sqrt(_x) + 
  _x*(NACA4_a1 + _x*(NACA4_a2 + _x*(NACA4_a3 + _x*NACA4_a4))) ;
  
  If ( _th > Pi ) _y = -_y ; EndIf
  
  _p = newp ;
  L1 = 1-_x ; L2 = _x ;
  NACA4_len = NACA4_len_le*L1*(2*L1-1) + NACA4_len_mp*4*L1*L2 +
  NACA4_len_te*L2*(2*L2-1) ;
  Point(_p) = {NACA4_ch*_x+NACA4_le_x,  
    _y*NACA4_ch*NACA4_th/20+NACA4_le_y, 
    NACA4_le_z, NACA4_len} ;
  NACA4_Points[_i] = _p ;
EndFor

For _i In {0:NACA4_nspl-1}
  _c = newc ;
  Spline(_c) = NACA4_Points[{_i*(NACA4_pps-1):(_i+1)*(NACA4_pps-1)}] ;
  NACA4_Splines[_i] = _c ;
EndFor

Return

Function NACA4_trailing_edge

// Generate the trailing edge of a symmetric NACA four series aerofoil
// composed of a set of splines, starting from a specified percentage
// of chord

//
// the following variables must be set on entry:
// NACA4_th     thickness in percent of chord
// NACA4_ch     aerofoil chord
// NACA4_le_x,y,z   leading edge coordinates
// NACA4_len_te length scale (trailing edge)
// NACA4_len_mc length scale (mid chord)
// NACA4_len_le length scale (leading edge)
// NACA4_nspl   number of splines on section
// NACA4_pps    number of points per spline
// NACA4_cut    point at which to cut the section, expressed as a 
//              fraction of chord

// The local scale length will be set using a quadratic which interpolates
// trailing edge, midpoint and leading edge scale lengths as a function of
// distance along the chord.

// On exit, the following variables will contain the details of
// the trailing edge
//
// NACA4_Points[] a list of the 2*NACA4_nspl*NACA4_pps points
//                on the section
// NACA4_Splines[] a list of the 2*NACA4_nspl splines
//
// these two lists are oriented so that they start at the trailing edge
// and move over the upper surface and around the lower surface to return
// to the trailing edge

NACA4_npts = NACA4_nspl*(NACA4_pps-1) ;

NACA4_Points[] = {} ;
NACA4_Splines[] = {} ;

_th0 = Acos(2*NACA4_cut - 1) ;

For _i In {0:NACA4_npts-1}
  _th = 2*Pi - _th0 + _i/(NACA4_npts-1)*(_th0) ;
  
  _x = 0.5*(1+Cos(_th)) ;
  _y = NACA4_a0*Sqrt(_x) + 
  _x*(NACA4_a1 + _x*(NACA4_a2 + _x*(NACA4_a3 + _x*NACA4_a4))) ;
  _y = -_y ;
  
  _p = newp ;
  L1 = 1-_x ; L2 = _x ;
  NACA4_len = NACA4_len_le*L1*(2*L1-1) + NACA4_len_mp*4*L1*L2 +
  NACA4_len_te*L2*(2*L2-1) ;
  Point(_p) = {NACA4_ch*_x+NACA4_le_x,  
    _y*NACA4_ch*NACA4_th/20+NACA4_le_y, 
    NACA4_le_z, NACA4_len} ;
  NACA4_Points[_i] = _p ;
EndFor

// Generate the trailing edge extension
_x = 1 ;
_y = NACA4_a0*Sqrt(_x) + 
    _x*(NACA4_a1 + _x*(NACA4_a2 + _x*(NACA4_a3 + _x*NACA4_a4))) ;
_dy = 0.5*NACA4_a0/Sqrt(_x) + 
    (NACA4_a1 + _x*(2*NACA4_a2 + _x*(3*NACA4_a3 + _x*4*NACA4_a4))) ;

_x = 1 - _y/_dy ; _y = 0 ;
_p = newp ;
Point(_p) = {NACA4_ch*_x+NACA4_le_x,  
	    _y*NACA4_ch*NACA4_th/20+NACA4_le_y, 
	    NACA4_le_z, NACA4_len_te} ;
NACA4_Points[NACA4_npts] = _p ;

For _i In {0:NACA4_npts-1}
  _th = _i/(NACA4_npts-1)*_th0 ;
  
  _x = 0.5*(1+Cos(_th)) ;
  _y = NACA4_a0*Sqrt(_x) + 
  _x*(NACA4_a1 + _x*(NACA4_a2 + _x*(NACA4_a3 + _x*NACA4_a4))) ;
  
  _p = newp ;
  L1 = 1-_x ; L2 = _x ;
  NACA4_len = NACA4_len_le*L1*(2*L1-1) + NACA4_len_mp*4*L1*L2 +
  NACA4_len_te*L2*(2*L2-1) ;
  Point(_p) = {NACA4_ch*_x+NACA4_le_x,  
    _y*NACA4_ch*NACA4_th/20+NACA4_le_y, 
    NACA4_le_z, NACA4_len} ;
  NACA4_Points[NACA4_npts+1+_i] = _p ;
EndFor

For _i In {0:2*NACA4_nspl-1}
  _c = newc ;
  Spline(_c) = NACA4_Points[{_i*(NACA4_pps-1):(_i+1)*(NACA4_pps-1)}] ;
  NACA4_Splines[_i] = _c ;
EndFor

NACA4_npts = 2*NACA4_npts+1 ;

Return

Function SectionJoin

// Join two aerofoil sections (from NACA4) to form a wing surface
//
// Variables which should be set on input:
//
// SJ_step     the number of points per spline on each section (must be
//             the same for both)
// SJ_points   the list of points forming the sections
// SJ_splines  the list of splines forming the sections
// SJ_offp1, SJ_offp2  the index of the trailing edge point for the first 
//             and second section in SJ_points
// SJ_offsp1, SJ_offsp2  the index of the first spline for the first 
//             and second section in SJ_points
// SJ_n_spline the number of splines per section

_c0 = newl ; _s0 = _c0 ;
Printf("Making trailing edge line (%g,%g)", 
       SJ_points[SJ_offp1],SJ_points[SJ_offp2]) ;
Line(_c0) = {SJ_points[SJ_offp1], SJ_points[SJ_offp2]} ;

_j = SJ_step-1 ;
For _i In {0:SJ_n_spline-2}
  _c1 = newc ;
  Printf("Making line (%g,%g)", SJ_points[SJ_offp1+_j],SJ_points[SJ_offp2+_j]) ;

  Line(_c1) = {SJ_points[SJ_offp1+_j],SJ_points[SJ_offp2+_j]} ;

  _c = newll ;
  Printf("Making loop(%g) (%g,%g,%g,%g)", _c,
    SJ_splines[SJ_offsp1+_i],_c1,-SJ_splines[SJ_offsp2+_i],-_c0) ;
  Line Loop(_c) = {SJ_splines[SJ_offsp1+_i],_c1,-SJ_splines[SJ_offsp2+_i],-_c0} ;
  s = news ;
  Ruled Surface(s) = {_c};
  _c0 = _c1 ; 
  _j = _j + SJ_step-1 ;
EndFor

_c = newll ;
Printf("Making loop(%g) = (%g,%g,%g,%g)", _c,
  SJ_splines[SJ_offsp1+SJ_n_spline-1],_s0,
  -SJ_splines[SJ_offsp2+SJ_n_spline-1],-_c0) ;
Line Loop(_c) = {SJ_splines[SJ_offsp1+SJ_n_spline-1],_s0,
  -SJ_splines[SJ_offsp2+SJ_n_spline-1],-_c0} ;
_s = news ;
Ruled Surface(_s) = {_c};

Return
