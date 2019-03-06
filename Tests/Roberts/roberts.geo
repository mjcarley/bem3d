Include "../GEO/naca.geo" ;

// Root chord
Root = 1.0 ;
// Wing span
Span = 4 ;
// Taper ratio, (Tip chord)/(Root chord)
Taper = 1/3 ;
// Leading edge sweep angle
Sweep = 41.81*Pi/180 ;

BaseLength = 2.0 ;

NACA4_nspl = 16 ;
NACA4_pps = 4 ;
NACA4_len_max = 0.2 ;
NACA4_len_min = 0.025 ;
NACA4_len_scale = 0.07 ;
NACA4_trunc = 1 ;

NACA4_th = 12 ;
NACA4_ch = Root*Taper ;
NACA4_le_x = Span/2*Sin(Sweep) ;
NACA4_le_y = 0.0 ;
NACA4_le_z = -Span/2 ;

Call NACA4 ;

ll = newll ;
Line Loop(ll) = NACA4_Splines[] ;
s = news ;
Plane Surface(s) = {-ll};

points[] = NACA4_Points[] ;
splines[] = NACA4_Splines[] ;

NACA4_th = 15 ;
NACA4_ch = Root ;
NACA4_le_x = 0.0 ;
NACA4_le_y = 0.0 ;
NACA4_le_z = 0.0 ;

Call NACA4 ;

points[] = {points[],NACA4_Points[]} ;
splines[] = {splines[],NACA4_Splines[]} ;

NACA4_ch = Root*Taper ;
NACA4_le_x = Span/2*Sin(Sweep) ;
NACA4_le_y = 0.0 ;
NACA4_le_z = Span/2 ;

Call NACA4 ;

ll = newll ;
Line Loop(ll) = NACA4_Splines[] ;
s = news ;
Plane Surface(s) = {ll};

points[] = {points[],NACA4_Points[]} ;
splines[] = {splines[],NACA4_Splines[]} ;

npts = NACA4_npts ;
nsplines = NACA4_nspl ;
SJ_step = NACA4_pps ;
SJ_points[] = points[] ;
SJ_splines[] = splines[] ;
SJ_n_spline = nsplines ;

SJ_offp1 = 0 ; SJ_offp2 = npts ;
SJ_offsp1 = 0 ; SJ_offsp2 = nsplines ;

Call SectionJoin ;

SJ_offp1 = npts ; SJ_offp2 = 2*npts ;
SJ_offsp1 = nsplines ; SJ_offsp2 = 2*nsplines ;

Call SectionJoin ;
