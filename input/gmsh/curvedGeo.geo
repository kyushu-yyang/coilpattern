// Gmsh project created on Sat Apr 23 15:08:04 2022
SetFactory("OpenCASCADE");

//+
Ellipse(1) = {1.89, 0, 0, 0.077, 0.055, 0.0, Pi*2.0};
Transfinite Curve {1} = 65 Using Progression 1;
//+
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/16} {
  Curve{1}; Layers{75}; 
}
//+
Physical Curve("boundary", 4) = {3, 1};
//+
Physical Surface("surface", 5) = {1};
