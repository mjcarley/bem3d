#include <stdio.h>

#include <math.h>

#include <glib.h>

#ifndef REAL_FLOAT
#ifndef REAL_DOUBLE
#define REAL_DOUBLE
#endif
#endif

#ifdef REAL_FLOAT
/* #warning "Compiling in single precision" */
#define REAL gfloat
#define REAL_ATAN2(x,y)  atan2f((x),(y))
#define REAL_ATAN(x)     atanf((x))
#define REAL_SIN(x)      sinf((x))
#define REAL_COS(x)      cosf((x))
#define REAL_LOG(x)      logf((x))
#define REAL_SQRT(x)     sqrtf((x))
#define NEWMAN_TRI                   newman_tri_f
#define NEWMAN_TRI_MP                newman_tri_mp_f
#define NEWMAN_TRI_SHAPE             newman_tri_shape_f
#define NEWMAN_QUAD                  newman_quad_f
#define NEWMAN_TRI_MOMENTS           newman_tri_moments_f
#define NEWMAN_TRI_GRADIENT          newman_tri_gradient_f
#define NEWMAN_TRI_SHAPE_GRADIENT    newman_tri_shape_gradient_f
#else
/* #warning "Compiling in double precision" */
#define REAL gdouble
#define REAL_ATAN2(x,y)  atan2((x),(y))
#define REAL_ATAN(x)     atan((x))
#define REAL_SIN(x)      sin((x))
#define REAL_COS(x)      cos((x))
#define REAL_LOG(x)      log((x))
#define REAL_SQRT(x)     sqrt((x))
#define NEWMAN_TRI         newman_tri
#define NEWMAN_TRI_MP      newman_tri_mp
#define NEWMAN_TRI_SHAPE   newman_tri_shape
#define NEWMAN_QUAD        newman_quad
#define NEWMAN_TRI_MOMENTS newman_tri_moments
#define NEWMAN_TRI_GRADIENT          newman_tri_gradient
#define NEWMAN_TRI_SHAPE_GRADIENT    newman_tri_shape_gradient
#endif /*REAL*/

const gdouble NEWMAN_FACTORIALS[] = {
  1.0000000000000000e+00, 1.0000000000000000e+00, 2.0000000000000000e+00,
  6.0000000000000000e+00, 2.4000000000000000e+01, 1.2000000000000000e+02,
  7.2000000000000000e+02, 5.0400000000000000e+03, 4.0320000000000000e+04,
  3.6288000000000000e+05, 3.6288000000000084e+06, 3.9916800000000030e+07,
  4.7900159999999928e+08, 6.2270207999999657e+09, 8.7178291200000168e+10,
  1.3076743679999983e+12, 2.0922789887999980e+13, 3.5568742809599812e+14,
  6.4023737057279940e+15, 1.2164510040883208e+17, 2.4329020081766400e+18,
  5.1090942171709784e+19, 1.1240007277776035e+21, 2.5852016738885062e+22,
  6.2044840173323914e+23, 1.5511210043331066e+25, 4.0329146112660538e+26,
  1.0888869450418268e+28, 3.0488834461171542e+29, 8.8417619937396265e+30,
  2.6525285981218941e+32, 8.2228386541778936e+33, 2.6313083693369503e+35,
  8.6833176188119946e+36, 2.9523279903960499e+38, 1.0333147966386149e+40,
  3.7199332678990102e+41, 1.3763753091226160e+43, 5.2302261746659953e+44,
  2.0397882081197180e+46, 8.1591528324790563e+47, 3.3452526613163927e+49,
  1.4050061177528771e+51, 6.0415263063373260e+52, 2.6582715747884193e+54,
  1.1962222086548034e+56, 5.5026221598119888e+57, 2.5862324151117192e+59,
  1.2413915592535963e+61, 6.0828186403425592e+62, 3.0414093201713019e+64,
  1.5511187532873646e+66, 8.0658175170944942e+67, 4.2748832840601228e+69,
  2.3084369733924540e+71, 1.2696403353658287e+73, 7.1099858780487453e+74,
  4.0526919504876951e+76, 2.3505613312829024e+78, 1.3868311854568818e+80,
  8.3209871127414438e+81, 5.0758021387722225e+83, 3.1469973260387477e+85,
  1.9826083154044198e+87, 1.2688693218588428e+89, 8.2476505920820880e+90,
  5.4434493907743010e+92, 3.6471110918189461e+94, 2.4800355424368993e+96,
  1.7112245242814200e+98, 1.1978571669969933e+100, 8.5047858856785845e+101,
  6.1234458376885870e+103, 4.4701154615127710e+105, 3.3078854415192728e+107,
  2.4809140811395329e+109, 1.8854947016660685e+111, 1.4518309202828782e+113,
  1.1324281178206039e+115, 8.9461821307828735e+116, 7.1569457046266771e+118,
  5.7971260207475536e+120, 4.7536433370131092e+122, 3.9455239697206508e+124,
  3.3142401345654752e+126, 2.8171041143804693e+128, 2.4227095383671897e+130,
  2.1077572983796612e+132, 1.8548264225740497e+134, 1.6507955160908248e+136,
  1.4857159644817212e+138, 1.3520015276784450e+140, 1.2438414054641828e+142,
  1.1567725070816412e+144, 1.0873661566567803e+146, 1.0329978488238504e+148,
  9.9167793487095964e+149, 9.6192759682489238e+151, 9.4268904488835978e+153,
  9.3326215443949336e+155, 9.3326215443942249e+157
} ;

/*binomials (n,k) for n=0,...,33*/

const gint NEWMAN_BINOMIALS[] = {
  1,
  1, 1,
  1, 2,  1,
  1, 3,  3,  1,
  1, 4,  6,  4,   1,
  1, 5, 10, 10,   5,   1,
  1, 6, 15, 20,  15,   6,   1,
  1, 7, 21, 35,  35,  21,   7,    1,
  1, 8, 28, 56,  70,  56,  28,    8,    1,
  1, 9, 36, 84, 126, 126,  84,   36,    9,    1,
  1,10, 45,120, 210, 252, 210,  120,   45,   10,   1,
  1,11, 55,165, 330, 462, 462,  330,  165,   55,  11,   1,
  1,12, 66,220, 495, 792, 924,  792,  495,  220,  66,  12,   1,
  1,13, 78,286, 715,1287,1716, 1716, 1287,  715, 286,  78,  13,  1,
  1,14, 91,364,1001,2002,3003, 3432, 3003, 2002,1001, 364,  91, 14,  1,
  1,15,105,455,1365,3003,5005, 6435, 6435, 5005,3003,1365, 455,105, 15, 1,
  1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1,
  1,17,136,680,2380,6188,12376,19448,24310,24310,19448,12376,6188,2380,680,136,17,1,
  1,18,153,816,3060,8568,18564,31824,43758,48620,43758,31824,18564,8568,3060,816,153,18,1,
  1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,3876,969,171,19,1,
  1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1,
  1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,352716,293930,203490,116280,54264,20349,5985,1330,210,21,1,
  1,22,231,1540,7315,26334,74613,170544,319770,497420,646646,705432,646646,497420,319770,170544,74613,26334,7315,1540,231,22,1,
  1,23,253,1771,8855,33649,100947,245157,490314,817190,1144066,1352078,1352078,1144066,817190,490314,245157,100947,33649,8855,1771,253,23,1,
  1,24,276,2024,10626,42504,134596,346104,735471,1307504,1961256,2496144,2704156,2496144,1961256,1307504,735471,346104,134596,42504,10626,2024,276,24,1,
  1,25,300,2300,12650,53130,177100,480700,1081575,2042975,3268760,4457400,5200300,5200300,4457400,3268760,2042975,1081575,480700,177100,53130,12650,2300,300,25,1,
  1,26,325,2600,14950,65780,230230,657800,1562275,3124550,5311735,7726160,9657700,10400600,9657700,7726160,5311735,3124550,1562275,657800,230230,65780,14950,2600,325,26,1,
  1,27,351,2925,17550,80730,296010,888030,2220075,4686825,8436285,13037895,17383860,20058300,20058300,17383860,13037895,8436285,4686825,2220075,888030,296010,80730,17550,2925,351,27,1,
  1,28,378,3276,20475,98280,376740,1184040,3108105,6906900,13123110,21474180,30421755,37442160,40116600,37442160,30421755,21474180,13123110,6906900,3108105,1184040,376740,98280,20475,3276,378,28,1,
  1,29,406,3654,23751,118755,475020,1560780,4292145,10015005,20030010,34597290,51895935,67863915,77558760,77558760,67863915,51895935,34597290,20030010,10015005,4292145,1560780,475020,118755,23751,3654,406,29,1,
  1,30,435,4060,27405,142506,593775,2035800,5852925,14307150,30045015,54627300,86493225,119759850,145422675,155117520,145422675,119759850,86493225,54627300,30045015,14307150,5852925,2035800,593775,142506,27405,4060,435,30,1,
  1,31,465,4495,31465,169911,736281,2629575,7888725,20160075,44352165,84672315,141120525,206253075,265182525,300540195,300540195,265182525,206253075,141120525,84672315,44352165,20160075,7888725,2629575,736281,169911,31465,4495,465,31,1,
  1,32,496,4960,35960,201376,906192,3365856,10518300,28048800,64512240,129024480,225792840,347373600,471435600,565722720,601080390,565722720,471435600,347373600,225792840,129024480,64512240,28048800,10518300,3365856,906192,201376,35960,4960,496,32,1,
  1,33,528,5456,40920,237336,1107568,4272048,13884156,38567100,92561040,193536720,354817320,573166440,818809200,1037158320,1166803110,1166803110,1037158320,818809200,573166440,354817320,193536720,92561040,38567100,13884156,4272048,1107568,237336,40920,5456,528,33,1
} ;

#define newman_binomial(m,k) (NEWMAN_BINOMIALS[(m)*((m)+1)/2+(k)])
#define newman_factorial(i) NEWMAN_FACTORIALS[(i)]
#define newman_r_block_start(h) ((h)*((h)+1)*((h)+2)/6)
#define newman_r_index_mnk(m,n,k) (((m)+(n))*((m)+(n)+1)/2 + (m))
#define newman_r_index(m,n,k) \
  (newman_r_block_start((m)+(n)+(k))+newman_r_index_mnk((m),(n),(k)))
#define newman_block_start(h) ((h)*((h)+1)/2)
#define newman_index_mn(m,n) (n)
#define newman_index(m,n) \
  (newman_block_start((m)+(n))+newman_index_mn((m),(n)))

#define _invert3x3(_Ai,_A)						\
  do {                                                                  \
    gdouble _det = 1.0/((_A)[0]*((_A)[8]*(_A)[4]-(_A)[7]*(_A)[5]) -     \
                        (_A)[3]*((_A)[8]*(_A)[1]-(_A)[7]*(_A)[2]) +     \
                        (_A)[6]*((_A)[5]*(_A)[1]-(_A)[4]*(_A)[2])) ;    \
    (_Ai)[0] =  _det*((_A)[8]*(_A)[4] - (_A)[7]*(_A)[5]) ;              \
    (_Ai)[1] = -_det*((_A)[8]*(_A)[1] - (_A)[7]*(_A)[2]) ;              \
    (_Ai)[2] =  _det*((_A)[5]*(_A)[1] - (_A)[4]*(_A)[2]) ;              \
                                                                        \
    (_Ai)[3] = -_det*((_A)[8]*(_A)[3] - (_A)[6]*(_A)[5]) ;              \
    (_Ai)[4] =  _det*((_A)[8]*(_A)[0] - (_A)[6]*(_A)[2]) ;              \
    (_Ai)[5] = -_det*((_A)[5]*(_A)[0] - (_A)[3]*(_A)[2]) ;              \
                                                                        \
    (_Ai)[6] =  _det*((_A)[7]*(_A)[3] - (_A)[6]*(_A)[4]) ;              \
    (_Ai)[7] = -_det*((_A)[7]*(_A)[0] - (_A)[6]*(_A)[1]) ;              \
    (_Ai)[8] =  _det*((_A)[4]*(_A)[0] - (_A)[3]*(_A)[1]) ;              \
  } while (0)


gint NEWMAN_TRI (REAL p[], REAL x1[], REAL x2[], REAL x3[],
		 REAL I[], REAL J[]) ;
gint NEWMAN_TRI_GRADIENT (REAL p[], REAL x1[], REAL x2[], REAL x3[],
			  REAL I[], REAL J[]) ;
gint NEWMAN_QUAD (REAL p[], 
		  REAL x1[], REAL x2[], 
		  REAL x3[], REAL x4[],
		  REAL I[], REAL J[]) ;
gint NEWMAN_TRI_MOMENTS (REAL x1[], REAL x2[],
			 REAL x3[], gint hmax, REAL Imn[]) ;

gint NEWMAN_TRI_SHAPE (REAL p[], REAL x1[], REAL x2[], REAL x3[],
		       REAL Imn[], gint hmax,
		       REAL I[], REAL J[]) ;

gint NEWMAN_TRI_SHAPE_GRADIENT (REAL p[], REAL x1[], REAL x2[], REAL x3[],
				REAL Imn[], gint hmax,
				REAL I[], REAL J[]) ;


/*
  Generate derivatives of the Laplace equation Green's function
  using the method of Tausch, J., `The fast multipole method for
  arbitrary Green's functions', Current Trends in Scientific
  Computing, 307--314, 2003.
 */

static void newman_laplace_tausch(REAL r[], gint H, REAL *dG)

{
  gint h, q, m, n, k, idx, jdx, kdx ;
  REAL DG[18000], R, R2 ;
  gint off[] = {0, 2600, 4900, 6924, 8695, 10235, 11565, 12705, 13674, 
		14490, 15170, 15730, 16185, 16549, 16835, 17055, 17220, 
		17340, 17424, 17480, 17515, 17535, 17545, 17549, 17550,
		17558} ;

  g_assert(H < 25) ;

  R2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ;
  R = sqrt(R2) ;

  for ( h = 0 ; h < 18000 ; h ++ ) DG[h] = 0.0 ;
  DG[0] = 1.0/R ; 
  for ( h = 0 ; h <= H+2 ; h ++ ) DG[off[h+1]] = -(2*h+1)*DG[off[h]]/R2 ;

  for ( q = H ; q >= 0 ; q -- ) {
    for ( h = 0 ; h <= H-q ; h ++ ) {
      for ( m = 0 ; m <= h ; m ++ ) {
	for ( n = 0 ; n <= h-m ; n ++ ) {
	  k = h-m-n ;
	  jdx = newman_r_index(m,n,k) ;

	  idx = newman_r_index(m+1,n,k) ; kdx = newman_r_index(m-1,n,k) ;
	  DG[off[q] + idx] =
	    m*DG[off[q+1]+kdx + 0] + r[0]*DG[off[q+1] + jdx] ;

	  idx = newman_r_index(m,n+1,k) ; kdx = newman_r_index(m,n-1,k) ;
	  DG[off[q] + idx] =
	    n*DG[off[q+1]+kdx + 0] + r[1]*DG[off[q+1] + jdx] ;

	  idx = newman_r_index(m,n,k+1) ; kdx = newman_r_index(m,n,k-1) ;
	  DG[off[q] + idx] =
	    k*DG[off[q+1]+kdx + 0] + r[2]*DG[off[q+1] + jdx] ;
	}	
      }
    }
  }

  for ( q = 0 ; q < newman_r_block_start(H+1) ; q ++ ) {
/*     if ( isnan(DG[q]) ) */
/*       g_error("%s: NaN at q==%d, H==%d, r=(%lg,%lg,%lg)", */
/* 	      __FUNCTION__, q, H, r[0], r[1], r[2]) ; */
    dG[q] = DG[q] ;
  }

  return ;
}

gint NEWMAN_TRI_MOMENTS (REAL x1[], REAL x2[],
			 REAL x3[], gint hmax, REAL Imn[])

{
  REAL z2pow[32], z3pow[32], n2pow[32], n3pow[32], A ;
  gint h, m, n, q, r ;

  for ( h = 0 ; h < newman_block_start(hmax+1) ; h ++ ) Imn[h] = 0.0 ;

  for ( (h = 1), (z2pow[0] = n2pow[0] = z3pow[0] = n3pow[0] = 1.0) ; 
	h <= hmax ; h ++ ) {
    z2pow[h] = z2pow[h-1]*x2[0] ; z3pow[h] = z3pow[h-1]*x3[0] ;
    n2pow[h] = n2pow[h-1]*x2[1] ; n3pow[h] = n3pow[h-1]*x3[1] ;
  }

  A = x2[0]*x3[1] - x2[1]*x3[0] ;

  for ( h = 0 ; h <= hmax ; h ++ ) {
    for ( m = 0 ; m <= h ; m ++ ) {
      n = h-m ;
      for ( q = 0 ; q <= m ; q ++ ) 
	for ( r = 0 ; r <= n ; r ++ ) 
	  Imn[newman_index(m,n)] += 
	    A*newman_binomial(m,q)*newman_binomial(n,r)*
	    z2pow[q]*z3pow[m-q]*n2pow[r]*n3pow[n-r]*
	    newman_factorial(q+r)*newman_factorial(h-q-r)/
	    newman_factorial(h+2) ;
    }
  }

  return 0 ;
}

/*
  Inputs: 
  p[]  position of field point in element-centred coordinates;

  x1, x2, x3: element vertices in element-centred coordinates (x1 must
  equal (0,0,0))

  Imn: element moments

  hmax: maximum order of moment, h=m+n

  Outputs:

  I, J: integrals over element
*/

static gint NEWMAN_TRI_MP (REAL p[], REAL x1[], REAL x2[],
			   REAL x3[], REAL Imn[], gint hmax,
			   REAL I[], REAL J[])

{
  REAL Rmn[512], sgn ;
  gint h, m, n ;
  REAL tol = 1e-6, err = 1e6, del ;

  g_assert( (x1[0] == x1[1]) && (x1[1] == x1[2]) && (x1[2] == 0.0) ) ;

  newman_laplace_tausch(p, hmax, Rmn) ;

  I[0] = I[1] = I[2] = J[0] = J[1] = J[2] = 0.0 ;
  for ( (h = 0), (sgn = 1.0) ; (h <= hmax-1) && (err > tol); 
	(h ++), (sgn = -sgn) ) {
    for ( m = 0 ; m <= h ; m ++ ) {
      n = h-m ;

      err = 0.0 ;
      J[0] -= (del = sgn*Imn[newman_index(m,n)]/
	       newman_factorial(m)/newman_factorial(n)*
	       Rmn[newman_r_index(m,n,0)]) ;
      err = MAX(err, fabs(del)) ;
      J[1] -= (del = sgn*Imn[newman_index(m+1,n)]/
	       newman_factorial(m)/newman_factorial(n)*
	       Rmn[newman_r_index(m,n,0)]) ;
      err = MAX(err, fabs(del)) ;
      J[2] -= (del = sgn*Imn[newman_index(m,n+1)]/
	       newman_factorial(m)/newman_factorial(n)*
	       Rmn[newman_r_index(m,n,0)]) ;
      err = MAX(err, fabs(del)) ;
      I[0] += (del = sgn*Imn[newman_index(m,n)]/
	       newman_factorial(m)/newman_factorial(n)*
	       Rmn[newman_r_index(m,n,1)]) ;
      err = MAX(err, fabs(del)) ;
      I[1] += (del = sgn*Imn[newman_index(m+1,n)]/
	       newman_factorial(m)/newman_factorial(n)*
	       Rmn[newman_r_index(m,n,1)]) ;
      I[2] += (del = sgn*Imn[newman_index(m,n+1)]/
	       newman_factorial(m)/newman_factorial(n)*
	       Rmn[newman_r_index(m,n,1)]) ;
      err = MAX(err, fabs(del)) ;
    }
  }

  return 0 ;
}

static gint newman_edge(REAL Rn, REAL Rnp1,
			REAL Sn, REAL Snp1,
			REAL *s, 
			REAL dn[], REAL dnp1[],
			REAL rn, REAL rnp1,
			REAL an, REAL anp1,
			REAL z,
			REAL I[], REAL J[])

{
  REAL Q, s1, s2, s3, c1, c2, c3, St, Ct, sn, u, U, tn ;

  tn = REAL_ATAN2(s[1], s[0]) ;
  St = REAL_SIN(tn) ; Ct = REAL_COS(tn) ;
  sn = REAL_SQRT(s[0]*s[0] + s[1]*s[1]) ;
  s1 = s[1]*Sn   - s[0]*dn[0]  *dn  [1] ; c1 = Rn  *z*s[0] ;
  s2 = s[1]*Snp1 - s[0]*dnp1[0]*dnp1[1] ; c2 = Rnp1*z*s[0] ;
  s3 = s1*c2 - s2*c1 ; c3 = c1*c2 + s1*s2 ;

  if ( c3 != 0.0 && z != 0.0 ) I[0] += REAL_ATAN2(s3,c3) ;

  if ( (Rn+Rnp1-sn) > 0.0 ) Q = REAL_LOG((Rn+Rnp1+sn)/(Rn+Rnp1-sn)) ;
  else Q = 0.0 ;

  J[0] += Q*(dn[0]*St - dn[1]*Ct) ;

  u = rn*REAL_COS(an - tn) ; U = rnp1*REAL_COS(anp1 - tn) ;

  I[1] += Q*St ; I[2] += Q*Ct ;

  Q = 0.5*(u*Rn - U*Rnp1 + (Rn*Rn - u*u)*Q) ;

  J[1] += Q*St ; J[2] += Q*Ct ;

  return 0 ;
}
			
gint NEWMAN_TRI (REAL p[], REAL x1[], REAL x2[], REAL x3[],
		 REAL I[], REAL J[])

{
  REAL R1, R2, R3, S1, S2, S3, r1, r2, r3, a1, a2, a3,
    d1[2], d2[2], d3[2] ;
  REAL s[2], z ;
  gint rt ;

  I[0] = I[1] = I[2] = J[0] = J[1] = J[2] = 0.0 ;
  z = p[2] ;

  r1 = (p[0]-x1[0])*(p[0]-x1[0]) +  (p[1]-x1[1])*(p[1]-x1[1]) ;
  r2 = (p[0]-x2[0])*(p[0]-x2[0]) +  (p[1]-x2[1])*(p[1]-x2[1]) ;
  r3 = (p[0]-x3[0])*(p[0]-x3[0]) +  (p[1]-x3[1])*(p[1]-x3[1]) ;

  R1 = REAL_SQRT(r1 + z*z) ; R2 = REAL_SQRT(r2 + z*z) ; 
  R3 = REAL_SQRT(r3 + z*z) ;
  S1 = (p[0]-x1[0])*(p[0]-x1[0]) + z*z ;
  S2 = (p[0]-x2[0])*(p[0]-x2[0]) + z*z ;
  S3 = (p[0]-x3[0])*(p[0]-x3[0]) + z*z ;

  r1 = REAL_SQRT(r1) ; r2 = REAL_SQRT(r2) ; r3 = REAL_SQRT(r3) ;

  a1 = REAL_ATAN2(p[1]-x1[1],p[0]-x1[0]) ;
  a2 = REAL_ATAN2(p[1]-x2[1],p[0]-x2[0]) ;
  a3 = REAL_ATAN2(p[1]-x3[1],p[0]-x3[0]) ;

  d1[0] = p[0] - x1[0] ; d1[1] = p[1] - x1[1] ;
  d2[0] = p[0] - x2[0] ; d2[1] = p[1] - x2[1] ;
  d3[0] = p[0] - x3[0] ; d3[1] = p[1] - x3[1] ;

  s[0] = x2[0]-x1[0] ; s[1] = x2[1]-x1[1] ;
  if ( (rt = newman_edge(R1, R2, S1, S2, s, d1, d2, r1, r2, a1, a2, z, I, J)) ) 
    return rt ;
  s[0] = x3[0]-x2[0] ; s[1] = x3[1]-x2[1] ;
  
  if ( (rt = newman_edge(R2, R3, S2, S3, s, d2, d3, r2, r3, a2, a3, z, I, J)) ) 
    return rt ;
  s[0] = x1[0]-x3[0] ; s[1] = x1[1]-x3[1] ;
  if ( (rt = newman_edge(R3, R1, S3, S1, s, d3, d1, r3, r1, a3, a1, z, I, J)) )
    return rt;

  I[1] = p[0]*I[0] + z*I[1] ;
  I[2] = p[1]*I[0] - z*I[2] ;

  J[0] -= z*I[0] ;
  J[1] = p[0]*J[0] - J[1] ;
  J[2] = p[1]*J[0] + J[2] ;

  return 0 ;
}

gint NEWMAN_QUAD (REAL p[], 
		  REAL x1[], REAL x2[], 
		  REAL x3[], REAL x4[],
		  REAL I[], REAL J[])

{
  REAL R1, R2, R3, R4, S1, S2, S3, S4, r1, r2, r3, r4,
    a1, a2, a3, a4, d1[2], d2[2], d3[2], d4[2] ;
  REAL s[2], z ;

  I[0] = I[1] = I[2] = J[0] = J[1] = J[2] = 0.0 ;
  z = p[2] ;

  r1 = (p[0]-x1[0])* (p[0]-x1[0]) +  (p[1]-x1[1])* (p[1]-x1[1]) ;
  r2 = (p[0]-x2[0])* (p[0]-x2[0]) +  (p[1]-x2[1])* (p[1]-x2[1]) ;
  r3 = (p[0]-x3[0])* (p[0]-x3[0]) +  (p[1]-x3[1])* (p[1]-x3[1]) ;
  r4 = (p[0]-x4[0])* (p[0]-x4[0]) +  (p[1]-x4[1])* (p[1]-x4[1]) ;

  R1 = REAL_SQRT(r1 + z*z) ; R2 = REAL_SQRT(r2 + z*z) ; 
  R3 = REAL_SQRT(r3 + z*z) ; R4 = REAL_SQRT(r4 + z*z) ;
  S1 = (p[0]-x1[0])*(p[0]-x1[0]) + z*z ;
  S2 = (p[0]-x2[0])*(p[0]-x2[0]) + z*z ;
  S3 = (p[0]-x3[0])*(p[0]-x3[0]) + z*z ;
  S4 = (p[0]-x4[0])*(p[0]-x4[0]) + z*z ;

  r1 = REAL_SQRT(r1) ; r2 = REAL_SQRT(r2) ; 
  r3 = REAL_SQRT(r3) ; r4 = REAL_SQRT(r4) ;

  a1 = REAL_ATAN2(p[1]-x1[1],p[0]-x1[0]) ;
  a2 = REAL_ATAN2(p[1]-x2[1],p[0]-x2[0]) ;
  a3 = REAL_ATAN2(p[1]-x3[1],p[0]-x3[0]) ;
  a4 = REAL_ATAN2(p[1]-x4[1],p[0]-x4[0]) ;

  d1[0] = p[0] - x1[0] ; d1[1] = p[1] - x1[1] ;
  d2[0] = p[0] - x2[0] ; d2[1] = p[1] - x2[1] ;
  d3[0] = p[0] - x3[0] ; d3[1] = p[1] - x3[1] ;
  d4[0] = p[0] - x4[0] ; d4[1] = p[1] - x4[1] ;
  
  s[0] = x2[0]-x1[0] ; s[1] = x2[1]-x1[1] ;
  newman_edge(R1, R2, S1, S2, s, d1, d2, r1, r2, a1, a2, z, I, J) ;
  s[0] = x3[0]-x2[0] ; s[1] = x3[1]-x2[1] ;
  newman_edge(R2, R3, S2, S3, s, d2, d3, r2, r3, a2, a3, z, I, J) ;
  s[0] = x4[0]-x3[0] ; s[1] = x4[1]-x3[1] ;
  newman_edge(R3, R4, S3, S4, s, d3, d4, r3, r4, a3, a4, z, I, J) ;
  s[0] = x1[0]-x4[0] ; s[1] = x1[1]-x4[1] ;
  newman_edge(R4, R1, S4, S1, s, d4, d1, r4, r1, a4, a1, z, I, J) ;

  I[1] = p[0]*I[0] + z*I[1] ;
  I[2] = p[1]*I[0] - z*I[2] ;

  J[0] -= z*I[0] ;
  J[1] = p[0]*J[0] - J[1] ;
  J[2] = p[1]*J[0] + J[2] ;

  return 0 ;
}

gint NEWMAN_TRI_SHAPE (REAL *p, REAL *x1, REAL *x2, REAL *x3,
		       REAL *Imn, gint hmax,
		       REAL *G, REAL *dG)

{
  REAL g[3], dg[3], A[9], Ai[9] ;

  if ( (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) >
       4*(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]) &&
       (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) >
       4*(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]) &&
       (Imn != NULL) ) {
    if ( NEWMAN_TRI_MP (p, x1, x2, x3, Imn, hmax, dg, g) ) return 1 ;
  } else {
    if ( NEWMAN_TRI (p, x1, x2, x3, dg, g) ) return 1 ;
  }

  A[0] = A[1] = A[2] = 1.0 ;
  A[3] = x1[0] ; A[4] = x2[0] ; A[5] = x3[0] ; 
  A[6] = x1[1] ; A[7] = x2[1] ; A[8] = x3[1] ; 

  _invert3x3(Ai, A) ;

  G[0] = Ai[0]*g[0] + Ai[1]*g[1] + Ai[2]*g[2] ;
  G[1] = Ai[3]*g[0] + Ai[4]*g[1] + Ai[5]*g[2] ;
  G[2] = Ai[6]*g[0] + Ai[7]*g[1] + Ai[8]*g[2] ;

  dG[0] = Ai[0]*dg[0] + Ai[1]*dg[1] + Ai[2]*dg[2] ;
  dG[1] = Ai[3]*dg[0] + Ai[4]*dg[1] + Ai[5]*dg[2] ;
  dG[2] = Ai[6]*dg[0] + Ai[7]*dg[1] + Ai[8]*dg[2] ;
  
  return 0 ;
}

static gint newman_edge_gradient(REAL Rn[], REAL Rnp1[],
				 REAL Sn[], REAL Snp1[],
				 REAL s[], 
				 REAL dn[], REAL dnp1[],
				 REAL rn[], REAL rnp1[],
				 REAL an[], REAL anp1[],
				 REAL z,
				 REAL I[], REAL J[])

{
  REAL P[4], Q[4], s1[4], s2[4], s3[4], c1[4], c2[4], c3[4], 
    St, Ct, sn, u[4], U[4], tn ;

  tn = REAL_ATAN2(s[1], s[0]) ;
  St = REAL_SIN(tn) ; Ct = REAL_COS(tn) ;
  sn = REAL_SQRT(s[0]*s[0] + s[1]*s[1]) ;

  s1[0] = s[1]*Sn[0]   - s[0]*dn[0]  *dn  [1] ; c1[0] = Rn[0]  *z*s[0] ;
  s2[0] = s[1]*Snp1[0] - s[0]*dnp1[0]*dnp1[1] ; c2[0] = Rnp1[0]*z*s[0] ;

  s1[1] = s[1]*Sn[1]   - s[0]*dn  [1] ; c1[1] = Rn[1]  *z*s[0] ;
  s2[1] = s[1]*Snp1[1] - s[0]*dnp1[1] ; c2[1] = Rnp1[1]*z*s[0] ;

  s1[2] = s[1]*Sn[2]   - s[0]*dn[0] ;   c1[2] = Rn[2]  *z*s[0] ;
  s2[2] = s[1]*Snp1[2] - s[0]*dnp1[0] ; c2[2] = Rnp1[2]*z*s[0] ;

  s1[3] = s[1]*Sn[3] ;
  c1[3] = Rn[3]  *z*s[0] + Rn[0]*s[0] ;
  s2[3] = s[1]*Snp1[3] ; 
  c2[3] = Rnp1[3]*z*s[0] + Rnp1[0]*s[0] ;

  s3[0] = s1[0]*c2[0] - s2[0]*c1[0] ; 
  c3[0] = c1[0]*c2[0] + s1[0]*s2[0] ;
  s3[1] = s1[0]*c2[1] + s1[1]*c2[0] - s2[0]*c1[1] - s2[1]*c1[0] ; 
  c3[1] = c1[0]*c2[1] + c1[1]*c2[0] + s1[0]*s2[1] + s1[1]*s2[0] ;
  s3[2] = s1[0]*c2[2] + s1[2]*c2[0] - s2[0]*c1[2] - s2[2]*c1[0] ; 
  c3[2] = c1[0]*c2[2] + c1[2]*c2[0] + s1[0]*s2[2] + s1[2]*s2[0] ;

  s3[3] = s1[0]*c2[3] + s1[3]*c2[0] - s2[0]*c1[3] - s2[3]*c1[0] ; 
  c3[3] = c1[0]*c2[3] + c1[3]*c2[0] + s1[0]*s2[3] + s1[3]*s2[0] ;

  if ( c3[0] != 0.0 && z != 0.0 ) {
    I[0] += REAL_ATAN2(s3[0],c3[0]) ;
    I[3] += (c3[0]*s3[1] - c3[1]*s3[0])/(s3[0]*s3[0] + c3[0]*c3[0]) ;
    I[6] += (c3[0]*s3[2] - c3[2]*s3[0])/(s3[0]*s3[0] + c3[0]*c3[0]) ;
    I[9] += (c3[0]*s3[3] - c3[3]*s3[0])/(s3[0]*s3[0] + c3[0]*c3[0]) ;
  }

  if ( (Rn[0]+Rnp1[0]-sn) > 0.0 ) {
    Q[0] = REAL_LOG((Rn[0]+Rnp1[0]+sn)/(Rn[0]+Rnp1[0]-sn)) ;
    Q[1] = -2.0*sn*(Rn[1] + Rnp1[1])/(Rn[0]+Rnp1[0]-sn)/(Rn[0]+Rnp1[0]+sn) ;
    Q[2] = -2.0*sn*(Rn[2] + Rnp1[2])/(Rn[0]+Rnp1[0]-sn)/(Rn[0]+Rnp1[0]+sn) ;
    Q[3] = -2.0*sn*(Rn[3] + Rnp1[3])/(Rn[0]+Rnp1[0]-sn)/(Rn[0]+Rnp1[0]+sn) ;
  } else Q[0] = Q[1] = Q[2] = Q[3] = 0.0 ;

  J[0] +=  Q[0]*(dn[0]*St - dn[1]*Ct) ;
  J[3] +=  Q[0]*St + (dn[0]*St - dn[1]*Ct)*Q[1] ;
  J[6] += -Q[0]*Ct + (dn[0]*St - dn[1]*Ct)*Q[2] ;
/*   J[9] += (dn[0]*St - dn[1]*Ct)*Q[3] ; */

  I[1]  += Q[0]*St ; I[2]  += Q[0]*Ct ;
  I[4]  += Q[1]*St ; I[5]  += Q[1]*Ct ;
  I[7]  += Q[2]*St ; I[8]  += Q[2]*Ct ;
  I[10] += Q[3]*St ; I[11] += Q[3]*Ct ;

  u[0] = rn[0]*REAL_COS(an[0] - tn) ; 
  u[1] = rn[1]*REAL_COS(an[0] - tn) - rn[0]*REAL_SIN(an[0]-tn)*an[1] ;
  u[2] = rn[2]*REAL_COS(an[0] - tn) - rn[0]*REAL_SIN(an[0]-tn)*an[2] ;
  u[3] = 0.0 ;

  U[0] = rnp1[0]*REAL_COS(anp1[0] - tn) ;
  U[1] = rnp1[1]*REAL_COS(anp1[0] - tn) - rnp1[0]*REAL_SIN(anp1[0]-tn)*anp1[1] ;
  U[2] = rnp1[2]*REAL_COS(anp1[0] - tn) - rnp1[0]*REAL_SIN(anp1[0]-tn)*anp1[2] ;
  U[3] = 0.0 ;

  P[0] = 0.5*(u[0]*Rn[0] - U[0]*Rnp1[0] + (Rn[0]*Rn[0] - u[0]*u[0])*Q[0]) ;
  P[1] = 0.5*(u[0]*Rn[1] + u[1]*Rn[0] - 
	      U[0]*Rnp1[1] - U[1]*Rnp1[0] + 
	      (Rn[0]*Rn[0] - u[0]*u[0])*Q[1] +
	      2.0*(Rn[0]*Rn[1] - u[0]*u[1])*Q[0]) ;
  P[2] = 0.5*(u[0]*Rn[2] + u[2]*Rn[0] - 
	      U[0]*Rnp1[2] - U[2]*Rnp1[0] + 
	      (Rn[0]*Rn[0] - u[0]*u[0])*Q[2] +
	      2.0*(Rn[0]*Rn[2] - u[0]*u[2])*Q[0]) ;
  P[3] = 0.5*(u[0]*Rn[3] + u[3]*Rn[0] - 
	      U[0]*Rnp1[3] - U[3]*Rnp1[0] + 
	      (Rn[0]*Rn[0] - u[0]*u[0])*Q[3] +
	      2.0*(Rn[0]*Rn[3] - u[0]*u[3])*Q[0]) ;

  J[1]  += P[0]*St ; J[2]  += P[0]*Ct ;
  J[4]  += P[1]*St ; J[5]  += P[1]*Ct ;
  J[7]  += P[2]*St ; J[8]  += P[2]*Ct ;
  J[10] += P[3]*St ; J[11] += P[3]*Ct ;

  return 0 ;
}

/*
  Compute the gradient of the three panel integrals on a triangle.

  Inputs:

  p:        field point in source coordinate system
  x1,x2,x3: nodes of triangular panel

  Outputs:

  I[]:      \int [1, x, y] \partial  (1/R)\partial z d A
  J[]:      \int [1, x, y] (1/R) d A

  Outputs are given as the basic integrals followed by the derivatives
  with respect to x, y and z.

*/
			
gint NEWMAN_TRI_GRADIENT (REAL p[], REAL x1[], REAL x2[], REAL x3[],
			  REAL I[], REAL J[])

{
  REAL R1[4], R2[4], R3[4], S1[4], S2[4], S3[4], 
    r1[4], r2[4], r3[4], a1[4], a2[4], a3[4], d1[2], d2[2], d3[2] ;
  REAL s[2], z ;
  gint rt ;

  I[0] = I[1] = I[2] = I[3] = I[4] = I[5] = I[6] = I[7] = I[8] = 
    I[9] = I[10] = I[11] = 
    J[0] = J[1] = J[2] = J[3] = J[4] = J[5] = J[6] = J[7] = J[8] = 
    J[9] = J[10] = J[11] = 0.0 ;
  z = p[2] ;

  r1[0] = (p[0]-x1[0])*(p[0]-x1[0]) +  (p[1]-x1[1])*(p[1]-x1[1]) ;
  r2[0] = (p[0]-x2[0])*(p[0]-x2[0]) +  (p[1]-x2[1])*(p[1]-x2[1]) ;
  r3[0] = (p[0]-x3[0])*(p[0]-x3[0]) +  (p[1]-x3[1])*(p[1]-x3[1]) ;

  R1[0] = REAL_SQRT(r1[0] + z*z) ;
  R1[1] = (p[0]-x1[0])/R1[0] ; R1[2] = (p[1]-x1[1])/R1[0] ; R1[3] = z/R1[0] ;

  R2[0] = REAL_SQRT(r2[0] + z*z) ;
  R2[1] = (p[0]-x2[0])/R2[0] ; R2[2] = (p[1]-x2[1])/R2[0] ; R2[3] = z/R2[0] ;

  R3[0] = REAL_SQRT(r3[0] + z*z) ;
  R3[1] = (p[0]-x3[0])/R3[0] ; R3[2] = (p[1]-x3[1])/R3[0] ; R3[3] = z/R3[0] ;

  S1[0] = (p[0]-x1[0])*(p[0]-x1[0]) + z*z ;
  S1[1] = 2*(p[0]-x1[0]) ; S1[2] = 0.0 ; S1[3] = 2*z ;
  
  S2[0] = (p[0]-x2[0])*(p[0]-x2[0]) + z*z ;
  S2[1] = 2*(p[0]-x2[0]) ; S2[2] = 0.0 ; S2[3] = 2*z ;

  S3[0] = (p[0]-x3[0])*(p[0]-x3[0]) + z*z ;
  S3[1] = 2*(p[0]-x3[0]) ; S3[2] = 0.0 ; S3[3] = 2*z ;

  r1[0] = REAL_SQRT(r1[0]) ; 
  r1[1] = (p[0]-x1[0])/r1[0] ; r1[2] = (p[1]-x1[1])/r1[0] ;

  r2[0] = REAL_SQRT(r2[0]) ; 
  r2[1] = (p[0]-x2[0])/r2[0] ; r2[2] = (p[1]-x2[1])/r2[0] ;

  r3[0] = REAL_SQRT(r3[0]) ;
  r3[1] = (p[0]-x3[0])/r3[0] ; r3[2] = (p[1]-x3[1])/r3[0] ;

  d1[0] = p[0] - x1[0] ; d1[1] = p[1] - x1[1] ;
  d2[0] = p[0] - x2[0] ; d2[1] = p[1] - x2[1] ;
  d3[0] = p[0] - x3[0] ; d3[1] = p[1] - x3[1] ;

  a1[0] = REAL_ATAN2(p[1]-x1[1],p[0]-x1[0]) ;
  a1[1] = -d1[1]/r1[0]/r1[0] ; a1[2] = d1[0]/r1[0]/r1[0] ;

  a2[0] = REAL_ATAN2(p[1]-x2[1],p[0]-x2[0]) ;
  a2[1] = -d2[1]/r2[0]/r2[0] ; a2[2] = d2[0]/r2[0]/r2[0] ;

  a3[0] = REAL_ATAN2(p[1]-x3[1],p[0]-x3[0]) ;
  a3[1] = -d3[1]/r3[0]/r3[0] ; a3[2] = d3[0]/r3[0]/r3[0] ;

  s[0] = x2[0]-x1[0] ; s[1] = x2[1]-x1[1] ;
  if ( (rt = newman_edge_gradient(R1, R2, S1, S2, s, d1, d2, r1, r2, 
				  a1, a2, z, I, J)) ) 
    return rt ;

  s[0] = x3[0]-x2[0] ; s[1] = x3[1]-x2[1] ;  
  if ( (rt = newman_edge_gradient(R2, R3, S2, S3, s, d2, d3, r2, r3, 
				  a2, a3, z, I, J)) ) 
    return rt ;

  s[0] = x1[0]-x3[0] ; s[1] = x1[1]-x3[1] ;
  if ( (rt = newman_edge_gradient(R3, R1, S3, S1, s, d3, d1, r3, r1, 
				  a3, a1, z, I, J)) )
    return rt;

  I[10] = p[0]*I[9] + I[1] + z*I[10] ;
  I[11] = p[1]*I[9] - I[2] - z*I[11] ;

  I[1]  =  p[0]*I[0] + z*I[1] ;
  I[2]  =  p[1]*I[0] - z*I[2] ;

  I[4]  =  I[0] + p[0]*I[3] + z*I[4] ;
  I[5]  =         p[1]*I[3] - z*I[5] ;

  I[7]  =         p[0]*I[6] + z*I[7] ;
  I[8]  =  I[0] + p[1]*I[6] - z*I[8] ;

  J[0] -= z*I[0] ;
  J[1]  =  p[0]*J[0] - J[1] ;
  J[2]  =  p[1]*J[0] + J[2] ;

  J[3] -= z*I[3] ;
  J[4]  =  J[0] + p[0]*J[3] - J[4] ;
  J[5]  =         p[1]*J[3] + J[5] ;

  J[6] -= z*I[6] ;
  J[7]  =         p[0]*J[6] - J[7] ;
  J[8]  =  J[0] + p[1]*J[6] + J[8] ;

  J[9]  = -I[0] ;
/*   J[9] -= z*I[9] + I[0] ; */
  J[10] = p[0]*J[9] - J[10] ;
  J[11] = p[1]*J[9] + J[11] ;

  return 0 ;
}

gint NEWMAN_TRI_SHAPE_GRADIENT (REAL *p, REAL *x1, REAL *x2, REAL *x3,
				REAL *Imn, gint hmax,
				REAL *G, REAL *dG)

{
  REAL g[12], dg[12], A[9], Ai[9] ;

  if ( NEWMAN_TRI_GRADIENT (p, x1, x2, x3, dg, g) ) return 1 ;

  A[0] = A[1] = A[2] = 1.0 ;
  A[3] = x1[0] ; A[4] = x2[0] ; A[5] = x3[0] ; 
  A[6] = x1[1] ; A[7] = x2[1] ; A[8] = x3[1] ; 

  _invert3x3(Ai, A) ;

  G[0] = Ai[0]*g[0] + Ai[1]*g[1] + Ai[2]*g[2] ;
  G[1] = Ai[3]*g[0] + Ai[4]*g[1] + Ai[5]*g[2] ;
  G[2] = Ai[6]*g[0] + Ai[7]*g[1] + Ai[8]*g[2] ;

  dG[0] = Ai[0]*dg[0] + Ai[1]*dg[1] + Ai[2]*dg[2] ;
  dG[1] = Ai[3]*dg[0] + Ai[4]*dg[1] + Ai[5]*dg[2] ;
  dG[2] = Ai[6]*dg[0] + Ai[7]*dg[1] + Ai[8]*dg[2] ;
  

  return 0 ;
}
