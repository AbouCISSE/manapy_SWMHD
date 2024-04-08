//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 2, 0, 1.0};
//+
Point(4) = {12, 2, 0, 1.0};
//+
Point(5) = {12, 0, 0, 1.0};
//+
Point(6) = {12, 0, 0, 1.0};
//+
Point(7) = {30, 0, 0, 1.0};
//+
Point(8) = {30, 20, 0, 1.0};
//+
Point(9) = {12, 20, 0, 1.0};
//+
Point(10) = {12, 10, 0, 1.0};
//+
Point(11) = {10, 10, 0, 1.0};
//+
Point(12) = {10, 20, 0, 1.0};
//+
Point(13) = {0, 20, 0, 1.0};
//+
Line(1) = {13, 12};
//+
Line(2) = {12, 11};
//+
Line(3) = {11, 10};
//+
Line(4) = {10, 9};
//+
Line(5) = {9, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 5};
//+
Line(8) = {5, 4};
//+
Line(9) = {4, 3};
//+
Line(10) = {3, 2};
//+
Line(11) = {2, 1};
//+
Line(12) = {1, 13};
//+
Physical Curve("3") = {1, 2, 3, 4, 5, 11, 9, 10, 8, 7};
//+
Physical Curve("1") = {12};
//+
Physical Curve("2") = {6};
//+
Curve Loop(1) = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1};
//+
Plane Surface(1) = {1};