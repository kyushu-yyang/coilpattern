// Gmsh project created on Sat Apr 23 15:08:04 2022
SetFactory("OpenCASCADE");

//+
Ellipse(1) = {1.89, 0, 0, 0.077, 0.055, 0, 2*Pi};
//+
Extrude {0, 0, 0.4} {
  Curve{1}; 
}
//+
Transfinite Curve {1, 3} = 55 Using Progression 1;
//+
Transfinite Curve {2} = 80 Using Progression 1;
//+
Physical Curve("bc", 4) = {3, 1};
//+
Physical Surface("surf", 5) = {1};
