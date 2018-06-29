ls = 0.4*lc;
r = 5; // Radius of the sphere
Br = r/2; // Cubic bulk dimension
Nc = 2*Br/lc; // Number of layers in the bulk edges
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//                   Structured internal cube 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Point(1) = {Br, Br, Br, ls};
Point(2) = {Br, -Br, Br, ls};
Point(3) = {Br, -Br, -Br, ls};
Point(4) = {Br, Br, -Br, ls};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Extrude {-2*Br, 0, 0} {
  Surface{1};Layers{Nc-1};
 }
 //
Transfinite Line {1, 2, 3, 4, 20, 11, 12, 16, 8, 9, 6, 7} =Nc Using Progression 1;
Transfinite Surface {13};
Transfinite Surface {1};
Transfinite Surface {25};
Transfinite Surface {26};
Transfinite Surface {17};
Transfinite Surface {21};
Surface Loop(29) = {25, 1, 13, 17, 21, 26};
Transfinite Volume{1} = {1, 2, 3, 4, 5, 6, 10, 14};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//                       Spherical boundary layer
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
p1 = newp; Point(p1)= {0.0,0.0,0.0,lc};
p2 = newp; Point(p2)= {r,0.0,0.0,ls};
p3 = newp; Point(p3)= {0,r,0.0,ls};
p4 = newp; Point(p4)= {-r,0,0.0,ls};
p5 = newp; Point(p5)= {0,-r,0.0,ls};
p6 = newp; Point(p6)= {0,0,-r,ls};
p7 = newp; Point(p7)= {0,0,r,ls};
c1 = newreg; Circle(c1) = {p2,p1,p3};
c2 = newreg; Circle(c2) = {p3,p1,p4};
c3 = newreg; Circle(c3) = {p4,p1,p5};
c4 = newreg; Circle(c4) = {p5,p1,p2};
c5 = newreg; Circle(c5) = {p3,p1,p6};
c6 = newreg; Circle(c6) = {p6,p1,p5};
c7 = newreg; Circle(c7) = {p5,p1,p7};
c8 = newreg; Circle(c8) = {p7,p1,p3};
c9 = newreg; Circle(c9) = {p2,p1,p7};
c10 = newreg; Circle(c10) = {p7,p1,p4};
c11 = newreg; Circle(c11) = {p4,p1,p6};
c12 = newreg; Circle(c12)= {p6,p1,p2};
Line Loop(13) = {c2,c8,-c10};
Line Loop(15) = {c10,c3,c7};
Line Loop(17) = {-c8,-c9,c1};
Line Loop(19) = {-c11,-c2,c5};
Line Loop(21) = {-c5,-c12,-c1};
Line Loop(23) = {-c3,c11,c6};
Line Loop(25) = {-c7,c4,c9};
Line Loop(27) = {-c4,c12,-c6};
Surface(28) = {13};
Surface(29) = {15};
Surface(30) = {17};
Surface(31) = {19};
Surface(32) = {21};
Surface(33) = {23};
Surface(34) = {25};
Surface(35) = {27};
Surface Loop(30) = {28, 31, 33, 29, 34, 35, 32, 30};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//                   Mesh volumes definition 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Volume(2) = {29,30};
Physical Volume(1) = {2,1};
Coherence Mesh;
Mesh.CharacteristicLengthFromPoints = 1;