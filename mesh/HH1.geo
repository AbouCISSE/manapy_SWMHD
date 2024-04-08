h_elem = 0.5;
Point(1) = {8, 8, 0, h_elem};
Point(2) = {-8, 8,   0, h_elem};
Point(3) = {-8, -8,   0, h_elem};
Point(4) = {8, -8, 0, h_elem};


Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1)={1};

Periodic Line {2} ={4};
Periodic Line {1} ={3};

Physical Line("5") = {4};
Physical Line("6") = {2};
Physical Line("7") = {3};
Physical Line("8") = {1};


Physical Surface("1") = {1} ;

