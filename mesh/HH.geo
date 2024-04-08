lc = 0.3;

Point(1) = {-8, -8, 0, lc};
Point(2) = {8, -8,  0, lc} ;
Point(3) = {8, 8, 0, lc} ;
Point(4) = {-8,  8, 0, lc} ;

//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Periodic Line {1} ={3};
Periodic Line {2} ={4};

Physical Line("1") = {4};
Physical Line("2") = {2};
Physical Line("3") = {1};
Physical Line("4") = {3};


Physical Surface("1") = {1} ;


