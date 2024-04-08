h_elem = .4;
Point(1) = {10, 10, 0, h_elem};
Point(2) = {0, 10,   0, h_elem};
Point(3) = {0, 0,   0, h_elem};
Point(4) = {10, 0., 0, h_elem};


Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Line("1") = {4};
Physical Line("2") = {2};
Physical Line("3") = {1};
Physical Line("4") = {3};


Physical Surface("1") = {1} ;




