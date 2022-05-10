// Gmsh project created on Sat Apr 23 15:08:04 2022
SetFactory("OpenCASCADE");

//+
Circle(1) = {0, 0, 0, 0.5, 0, Pi*2 };
//+
Transfinite Curve {1} = 60 Using Progression 1;
//+
Extrude {0, 0, 2} {
  Curve{1}; Layers {80}; 
}
//+
Physical Curve("boundary", 5) = {3, 1};
Physical Surface("surface", 4) = {1};
