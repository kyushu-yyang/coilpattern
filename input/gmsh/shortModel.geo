// Gmsh project created on Sat Apr 23 15:08:04 2022
SetFactory("OpenCASCADE");

// curvature length
curvature = Pi/16;

// 1st layer
Ellipse(1) = {1.89, 0, 0, 0.076135, 0.054135, 0, 2*Pi};

// 2nd layer
Ellipse(2) = {1.89, 0, 0, 0.077325, 0.055325, 0, 2*Pi};

// 3rd layer
Ellipse(3) = {1.89, 0, 0, 0.078515, 0.056515, 0, 2*Pi};

// 4th layer
Ellipse(4) = {1.89, 0, 0, 0.079705, 0.057705, 0, 2*Pi};

Extrude {{0, 1, 0}, {0, 0, 0}, curvature} { Curve{1}; }
Extrude {{0, 1, 0}, {0, 0, 0}, curvature} { Curve{2}; }
Extrude {{0, 1, 0}, {0, 0, 0}, curvature} { Curve{3}; }
Extrude {{0, 1, 0}, {0, 0, 0}, curvature} { Curve{4}; }

//+
Physical Curve("bc", 13) = {8, 6, 10, 12, 1, 2, 3, 4};
//+
Physical Surface("surface", 14) = {1, 3, 4, 2};
//+
Transfinite Curve {8, 6, 10, 12, 4, 3, 2, 1} = 25 Using Progression 1;
//+
Transfinite Curve {11, 9, 7, 5} = 25 Using Progression 1;
