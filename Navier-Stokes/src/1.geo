// Gmsh project created on Fri May 31 12:36:23 2019
//+
SetFactory("Built-in");
//+
SetFactory("OpenCASCADE");
//+
Point(5) = {-20, -3, 0, 1.0};
//+
Point(6) = {-20, 3, 0, 1.0};
//+
Point(7) = {20, 3, 0, 1.0};
//+
Point(8) = {20, -3, 0, 1.0};
//+
Point(9) = {3, -3, 0, 1.0};
//+
Point(10) = {3, 3, 0, 1.0};
//+
Point(11) = {-3, 3, 0, 1.0};
//+
Point(12) = {-3, -3, 0, 1.0};
//+
Circle(11) = {0, -3, 0, 3, Pi, 2*Pi};
//+
Circle(12) = {0, 3, 0, 3, 0, Pi};
//+
Circle(13) = {0, 0, 0, 1, 0, 2*Pi};
//+
Line(1) = {5,6} ;
Line(2) = {6, 11} ;
Line(3) = {10, 7} ;
Line(4) = {7, 8} ;
Line(5) = {8, 9} ;
Line(6) = {12, 5} ;
// Granice domene
Line Loop(1) = {1, 2, 12, 3, 4, 5, 11, 6};
Line Loop(2) = {13};

Plane Surface(1) = {1,2} ;

