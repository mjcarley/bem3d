Include "./dimensions.geo" ;

Point( 1) = {-ls, -ws, -th, lc} ;
Point( 2) = { ls, -ws, -th, lc} ;
Point( 3) = { ls,  ws, -th, lc} ;
Point( 4) = {-ls,  ws, -th, lc} ;

Point( 5) = {-ls, -ws,   0, lc} ;
Point( 6) = { ls, -ws,   0, lc} ;
Point( 7) = { ls,  ws,   0, lc} ;
Point( 8) = {-ls,  ws,   0, lc} ;

Point( 9) = {-lp, -wp,   0, lc} ;
Point(10) = { lp, -wp,   0, lc} ;
Point(11) = { lp,  wp,   0, lc} ;
Point(12) = {-lp,  wp,   0, lc} ;

Line( 1) = {  1,  2} ;
Line( 2) = {  2,  3} ;
Line( 3) = {  3,  4} ;
Line( 4) = {  4,  1} ;

Line( 5) = {  5,  6} ;
Line( 6) = {  6,  7} ;
Line( 7) = {  7,  8} ;
Line( 8) = {  8,  5} ;

Line( 9) = {  1,  5} ;
Line(10) = {  2,  6} ;
Line(11) = {  3,  7} ;
Line(12) = {  4,  8} ;

Line(13) = {  9, 10} ;
Line(14) = { 10, 11} ;
Line(15) = { 11, 12} ;
Line(16) = { 12,  9} ;


//+
Curve Loop(1) = {7, 8, 5, 6};
//+
Curve Loop(2) = {15, 16, 13, 14};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {5, -10, -1, 9};
//+
Plane Surface(2) = {-3};
//+
Curve Loop(4) = {10, 6, -11, -2};
//+
Plane Surface(3) = {-4};
//+
Curve Loop(5) = {11, 7, -12, -3};
//+
Plane Surface(4) = {-5};
//+
Curve Loop(6) = {8, -9, -4, 12};
//+
Plane Surface(5) = {-6};
//+
Curve Loop(7) = {1, 2, 3, 4};
//+
Plane Surface(6) = {-7};
