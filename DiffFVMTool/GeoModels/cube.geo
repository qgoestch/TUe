lc = lc;
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {10,0.0,0.0,lc};
Point(3) = {10,10,0.0,lc};
Point(4) = {0,10,0.0,lc};
Line(1) = {4,3};
Line(2) = {3,2};
Line(3) = {2,1};
Line(4) = {1,4};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
tmp[] = Extrude {0,0.0,10} {
  Surface{6};
};
Physical Volume(1) = tmp[1];
Mesh.CharacteristicLengthFromPoints = 1;




