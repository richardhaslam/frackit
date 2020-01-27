SetFactory("OpenCASCADE");

// Bottom Volume
Point(1) = {0, 0, 0, 1.0};
Point(2) = {50, 0, 25, 1.0};
Point(3) = {100, 0, 0, 1.0};
Point(4) = {150, 0, 0, 1.0};
Bezier(1) = {1, 2, 3, 4};

Extrude {0, 100, 0} { Curve{1}; }
Extrude {0, 0, 25} { Surface{1}; }

// Center Volume
Point(11) = {0, 0, 55, 1.0};
Point(12) = {50, 0, 55, 1.0};
Point(13) = {120, 0, 40, 1.0};
Point(14) = {150, 0, 55, 1.0};
Bezier(13) = {11, 12, 13, 14};

Extrude {0, 100, 0} { Curve{13}; }
Line(17) = {11, 7};
Line(18) = {14, 9};
Line(19) = {16, 10};
Line(20) = {15, 8};
Curve Loop(8) = {7, -20, -14, 17};
Plane Surface(8) = {8};
Curve Loop(9) = {17, 11, -18, -13};
Plane Surface(9) = {9};
Curve Loop(10) = {10, -19, -15, 18};
Plane Surface(10) = {10};
Curve Loop(12) = {19, -12, -20, 16};
Plane Surface(11) = {12};
Surface Loop(2) = {8, 11, 10, 7, 9, 6};
Volume(2) = {2};

// Top Volume
Point(18) = {0, 0, 75, 1.0};
Point(19) = {50, 0, 75, 1.0};
Point(20) = {120, 0, 60, 1.0};
Point(21) = {150, 0, 65, 1.0};
Bezier(21) = {18, 19, 20, 21};
Extrude {0, 100, 0} { Curve{21}; }
Line(25) = {18, 11};
Line(26) = {21, 14};
Line(27) = {23, 16};
Line(28) = {22, 15};
Curve Loop(16) = {28, -14, -25, 22};
Plane Surface(15) = {16};
Curve Loop(17) = {25, 13, -26, -21};
Plane Surface(16) = {17};
Curve Loop(18) = {26, 15, -27, -23};
Plane Surface(17) = {18};
Curve Loop(19) = {27, -16, -28, 24};
Plane Surface(18) = {19};
Surface Loop(3) = {14, 15, 18, 17, 16, 7};
Volume(3) = {3};
