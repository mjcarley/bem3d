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
Line Loop(1) = {6, 7, 8, 9};
//+
Plane Surface(1) = {-1};
//+
Recursive Delete {
  Curve{1}; Curve{18}; Curve{4}; Curve{19}; Curve{3}; Curve{21}; Curve{2}; Curve{20}; 
}
//+
Recursive Delete {
  Curve{12}; Curve{13}; Curve{10}; Curve{11}; 
}
//+
Recursive Delete {
  Point{5}; Point{6}; Point{7}; Point{8}; 
}
