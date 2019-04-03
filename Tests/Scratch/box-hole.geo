lc = 0.1 ;

dx = 1.0 ;
dy = 1.0 ;
dz = 0.125 ;

dxp = 0.25 ;
dyp = 0.125 ;

Point(1) = {-dx, -dy, -dz, lc} ;
Point(2) = { dx, -dy, -dz, lc} ;
Point(3) = { dx,  dy, -dz, lc} ;
Point(4) = {-dx,  dy, -dz, lc} ;

Point(5) = {-dx, -dy,  0, lc} ;
Point(6) = { dx, -dy,  0, lc} ;
Point(7) = { dx,  dy,  0, lc} ;
Point(8) = {-dx,  dy,  0, lc} ;

Point(9) =  {-dxp, -dyp,  0, lc} ;
Point(10) = { dxp, -dyp,  0, lc} ;
Point(11) = { dxp,  dyp,  0, lc} ;
Point(12) = {-dxp,  dyp,  0, lc} ;


//+
Line(1) = {8, 7};
//+
Line(2) = {7, 6};
//+
Line(3) = {6, 5};
//+
Line(4) = {5, 8};
//+
Line(5) = {12, 12};
//+
Line(6) = {12, 11};
//+
Line(7) = {11, 10};
//+
Line(8) = {10, 9};
//+
Line(9) = {9, 12};
//+
Line(10) = {1, 2};
//+
Line(11) = {2, 3};
//+
Line(12) = {3, 4};
//+
Line(13) = {4, 1};
//+
Line(14) = {5, 1};
//+
Line(15) = {8, 4};
//+
Line(16) = {7, 3};
//+
Line(17) = {6, 2};
//+
Line(18) = {8, 12};
//+
Line(19) = {9, 5};
//+
Line(20) = {11, 7};
//+
Line(21) = {10, 6};
//+
//Line Loop(1) = {6, 7, 8, 9};
//+
//Plane Surface(1) = {-1};
//+
Line Loop(2) = {1, -20, -6, -18};
//+
Plane Surface(2) = {-2};
//+
Line Loop(3) = {18, -9, 19, 4};
//+
Plane Surface(3) = {-3};
//+
Line Loop(4) = {19, -3, -21, 8};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {21, -2, -20, 7};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {12, 13, 10, 11};
//+
Plane Surface(6) = {-6};
//+
Line Loop(7) = {17, 11, -16, 2};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {12, -15, 1, 16};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {15, 13, -14, 4};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {3, 14, 10, -17};
//+
Plane Surface(10) = {10};
