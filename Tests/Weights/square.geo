lc = 0.03125 ;

dx = 0.3 ;
dy = 0.7 ;

Point(1) = {-dx, -dy, 0, lc} ;
Point(2) = { dx, -dy, 0, lc} ;
Point(3) = { dx,  dy, 0, lc} ;
Point(4) = {-dx,  dy, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;


//+
Line Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
