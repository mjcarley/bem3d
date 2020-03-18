Include "./dimensions.geo" ;

Point(1) = {-lp, -wp,   0, lc} ;
Point(2) = { lp, -wp,   0, lc} ;
Point(3) = { lp,  wp,   0, lc} ;
Point(4) = {-lp,  wp,   0, lc} ;

Line( 1) = {  1,  2} ;
Line( 2) = {  2,  3} ;
Line( 3) = {  3,  4} ;
Line( 4) = {  4,  1} ;

Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};

