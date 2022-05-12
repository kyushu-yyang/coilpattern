SetFactory("OpenCASCADE");

Point(1) = {1.89-0.08, -0.06, 0, 1.0};
Point(2) = {1.89+0.08, -0.06, 0, 1.0};
Point(3) = {1.89+0.08,  0.06, 0, 1.0};
Point(4) = {1.89-0.08,  0.06, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Curve {1} = 50 Using Progression 1;
Transfinite Curve {2} = 50 Using Progression 1;
Transfinite Curve {3} = 50 Using Progression 1;
Transfinite Curve {4} = 50 Using Progression 1;

Extrude {{0, 1, 0}, {0, 0, 0}, Pi/16} { Curve{1,2,3,4}; Layers{75}; }

Physical Curve("boundary", 13) = {11, 12, 7, 9, 3, 4, 1, 2};
Physical Surface("surface", 14) = {3, 4, 1, 2};

