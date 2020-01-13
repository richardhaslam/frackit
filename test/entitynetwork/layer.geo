SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0, 1.0};
Point(2) = {50, 0, 25, 1.0};
Point(3) = {100, 0, 0, 1.0};
Point(4) = {150, 0, 0, 1.0};
Bezier(1) = {1, 2, 3, 4};

Extrude {0, 100, 0} { Curve{1}; }
Extrude {0, 0, 25} { Surface{1}; }
