// Discretisation
h = 0.03;

g2elab.mipse.common.object.Point(1) = {2,0,1,h};
g2elab.mipse.common.object.Point(2) = {2,0,-1,h};
g2elab.mipse.common.object.Point(3) = {-2,0,-1,h};
g2elab.mipse.common.object.Point(4) = {-2,0,1,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {4,1,2,3};
Plane Surface(6) = {5};
